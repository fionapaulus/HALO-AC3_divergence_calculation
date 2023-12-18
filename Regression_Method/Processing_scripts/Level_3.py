import warnings

import datetime
import os
import glob
from tqdm import tqdm

import numpy as np
import xarray as xr
import pandas as pd
import netCDF4
from configparser import ConfigParser

import metpy.calc as mpcalc
import metpy.interpolate as mpinterp
from metpy.units import units

warnings.filterwarnings("ignore")

config = ConfigParser(inline_comment_prefixes="#")                                          # read config file
config.read("config.cfg")

alt_var    = config["LEVEL 3"]["altitude_varname"]                                          # get variable names from config file
lat_var    = config["LEVEL 3"]["latitude_varname"]
lon_var    = config["LEVEL 3"]["longitude_varname"]
time_var   = config["LEVEL 3"]["time_varname"]
p_var      = config["LEVEL 3"]["pressure_varname"]
rh_var     = config["LEVEL 3"]["relative_humidity_varname"]
ta_var     = config["LEVEL 3"]["air_temperature_varname"]
dp_var     = config["LEVEL 3"]["dew_point_varname"]
q_var      = config["LEVEL 3"]["specific_humidity_varname"]
qv_var     = config["LEVEL 3"]["water_vapor_q_varname"]
qsat_var   = config["LEVEL 3"]["saturation_q_varname"]

def interp_along_height(dataset, height_limit=10000, vertical_spacing=10, max_gap=50, interp_method="bin"):
    """
    Function to interpolate all values along the gpsalt dimension of a netCDF dataset
    to a specified vertical spacing (10 m default) upto a given height level (10 km default)
    Given dataset must have data variables along the gpsalt dimension

    Arguments:    
        dataset          (xr.DataSet): Dataset with variables along "gpsalt" dimension
        height_limit     (float):      Altitude up to which interpolation is carried 
                                       out (default=10km)
        vertical_spacing (float):      Vertical spacing to which values are interpolated 
                                       (default=10m)
        max_gap          (float):      No interpolation if gap between two datapoints 
                                       is > max_gap (default=50m)
        interp_method    (string):     Interpolation method ("bin" or "linear_interpolate")

    Returns:
        new_interpolated_ds (xr.DataSet): New dataset with given datasets variables 
                                          interpolated at given vertical_spacing up to 
                                          given height_limit
    """

    interpolation_grid = np.arange(0, height_limit+vertical_spacing, vertical_spacing)      # in meters
    interpolation_bins = np.arange(-2, height_limit+vertical_spacing+2,                     # Bins len(interpolation_grid)+1;
                                   vertical_spacing).astype("int")                          # (a,b]; (meters)

    if interp_method == "linear_interpolate":                                               # linear interpolation method
        
        new_index           = np.arange(0, height_limit+vertical_spacing, vertical_spacing)
        new_interpolated_ds = dataset.interp(gpsalt=interpolation_grid)

    elif interp_method == "bin":                                                            # bin interpolation method
        
        new_interpolated_ds = dataset.astype(np.float32).groupby_bins(                      # groupby does not bin lat,lon and 
                              alt_var, interpolation_bins,                                  # time since they are coordinates
                              labels=interpolation_grid,
                              restore_coord_dims=True).mean()
        dataset[time_var]    = ([alt_var], dataset[time_var].values.astype(float))
        
        for coords in [lat_var, lon_var, time_var]:                                         # adding them as extra variables
            new_interpolated_ds[coords] = (dataset[coords].groupby_bins(alt_var,
                                           interpolation_bins, labels=interpolation_grid,
                                           restore_coord_dims=False).mean())

        new_interpolated_ds["sonde_id"] = dataset["sonde_id"]

        new_interpolated_ds = new_interpolated_ds.transpose()
        new_interpolated_ds = new_interpolated_ds.rename({alt_var+"_bins": alt_var})
        new_interpolated_ds = new_interpolated_ds.interpolate_na(
                              alt_var, max_gap=max_gap, use_coordinate=True)

        new_interpolated_ds[time_var]          = ([alt_var], pd.DatetimeIndex(
                                                  new_interpolated_ds.time.values))
        new_interpolated_ds.encoding[time_var] = {"units": "seconds since 2020-01-01",
                                                  "dtype": "int32", "_FillValue": np.iinfo("int32").max}
        new_interpolated_ds                    = new_interpolated_ds.rename(
                                                 {time_var: "interpolated_time"})

    return new_interpolated_ds

def interpolate_for_grid(file_path, height_limit=10000, vertical_spacing=10):
    """
    Function to interpolate a dataset with Level-2 data, in the format for Level-3 
    gridding
    
    Arguments:
        file_path           (string):  Path to Level-2 .nc file
        height_limit        (integer): Height up to which the interpolation is performed
                                       (default=10km)
        vertical_spacing    (integer): Vertical grid spacing (default 10m)
    
    Returns:
        interpolated_dataset (xr.DataSet): interpolated dataset
    """

    dataset = xr.open_dataset(file_path).swap_dims({time_var: alt_var})                     # creating a Dataset of necessary 
                                                                                            # shape to interpolate

    if q_var not in list(dataset.data_vars):                                                # if specific humidity
        w_s = mpcalc.saturation_mixing_ratio(dataset[p_var].values * units.Pa,              # is not already included in the 
                     dataset[ta_var].values * units.kelvin).magnitude                       # dataset, it is added here
        w   = dataset[rh_var].values * w_s
        q   = w / (1 + w)
        dataset[q_var] = (dataset[p_var].dims, q)
        dataset[q_var].attrs["units"] = "kg/kg"
        dataset[q_var].attrs["long_name"] = "specific humidity"

    if qsat_var not in list(dataset.data_vars):                                             # adding saturation specific humidity
        esl   = float(config["CONSTANTS"]["es0"]) * np.exp(float(config["CONSTANTS"]["at"]) # Tetens formula for saturation vapor
                * (dataset[dp_var].values - float(config["CONSTANTS"]["tmelt"]))            # pressure over water, in Pa
                / (dataset[dp_var].values - float(config["CONSTANTS"]["bt"])))   
        esi   = float(config["CONSTANTS"]["es0"])*np.exp(float(config["CONSTANTS"]["at_i"]) # Tetens formula for saturation vapor 
                * (dataset[dp_var].values - float(config["CONSTANTS"]["tmelt"]))            # pressure over ice, in Pa
                / (dataset[dp_var].values - float(config["CONSTANTS"]["bt_i"])))   

        fac_phase = (dataset[ta_var].values - float(config["CONSTANTS"]["tdn"])) / (        # calculation of the phase 
                    float(config["CONSTANTS"]["tup"]) - float(config["CONSTANTS"]["tdn"]))
        fac_phase = np.where( fac_phase<0., 0., fac_phase)
        fac_phase = np.where( fac_phase>1., 1., fac_phase)

        qsat = (float(config["CONSTANTS"]["rd"]) / float(config["CONSTANTS"]["rv"])         # calculation of saturation 
               ) * esl / (dataset[p_var].values - (1-float(config["CONSTANTS"]["rd"])       # specific humidity
                / float(config["CONSTANTS"]["rv"])) * esl)
        dataset[qsat_var] = (dataset[p_var].dims, qsat)
        dataset.qsat.attrs["units"] = "kg/kg"
        dataset.qsat.attrs["long_name"] = "saturation specific humidity"

    if qv_var not in list(dataset.data_vars):                                               # adding water vapor specific humidity
        qv = np.where(dataset.qsat.values>=0., dataset.qsat.values, None)                   
        dataset[qv_var] = (dataset[p_var].dims, np.float32(qv))
        dataset.qv.attrs["units"] = "kg/kg"
        dataset.qv.attrs["long_name"] = "water vapor specific humidity"

    interpolated_dataset = interp_along_height(dataset, height_limit=height_limit,          # Interpolation of all values along 
                                vertical_spacing=vertical_spacing,                          # the height dimension
                                max_gap=float(config["GRID SETTINGS"]["max_gap_fill"]),
                                interp_method=config["GRID SETTINGS"]["method"])

    global_attrs = {                                                                        # set global attributes for output
        "title": "Dropsondes processing Level 3",                                           # files
        "platform_id":           Platform,
        "instrument_id":         config["INSTRUMENT"]["instrument_id"],
        "product_id":            "interpolated",
        "author":                config["AUTHOR"]["name"],
        "author_email":          config["AUTHOR"]["email"],
        "author_institute":      config["AUTHOR"]["institute"],
        "creation_time":         str(datetime.datetime.utcnow()) + " UTC"
    }
    interpolated_dataset.attrs = global_attrs                                               # add global attributes

    for var in dataset.variables:                                                           # add attributes to all variables
        if var != time_var:
            interpolated_dataset[var].attrs = dataset[var].attrs
    
    return interpolated_dataset

def check_for_faults(list_of_files, height_limit=10000, vertical_spacing=10):
    """
    Function to check if the recorded gps Altitude matches the given height range or if 
    the recording failed in this regard.

    Arguments:
        list_of_files: List of filepaths for all recorded dropsondes of the flight
        height_limit: Height up to which the interpolation is performed (default 10km)
        vertical spacing: Vertical grid spacing (default 10m)

    Returns:
        output_list_of_files (list of strings): list of all filepaths for all recorded 
                                                dropsondes with reasonable gps Altitudes

    Records:
        output_log (.txt): file in Gridded-folder recording the filepaths of the
                           discarded sondes
    """

    output_list_of_files = []

    if os.path.exists(archive + "Level_3/processing_log.txt"):                              # check if the ouput file already
        os.remove(archive + "Level_3/processing_log.txt")                                   # exists and remove old version

    output_log = open(archive + "Level_3/processing_log.txt", "w")                          # create new output file
    output_log.write(f"Processing Log for gridding of dropsondes data in\n")                # add general information
    output_log.write(f"{archive}\n\n")                                                      
    output_log.write(f"created by:\n")
    output_log.write("  " + config["AUTHOR"]["name"] + "\n")
    output_log.write("  " + config["AUTHOR"]["institute"] + "\n")
    output_log.write("  " + config["AUTHOR"]["email"] + "\n\n")
    output_log.write(f"{datetime.datetime.now()}\n\n")
    
    for idx in np.arange(len(list_of_files)):
        dataset = xr.open_dataset(list_of_files[idx])

        if np.sum((dataset[alt_var]>-5)&(dataset[alt_var]<(height_limit+vertical_spacing+5)))>0:
            output_list_of_files = np.append(output_list_of_files, list_of_files[idx])
        else:
            output_log.write(f"altitude stored in {list_of_files[idx][-21:]} does not" +
                               "include data betweeen the ground and the defined height" +
                               "limit and was discarded for gridding\n\n")
            print(f"altitude stored stored in {list_of_files[idx][-21:]} does not " +
                   "include data between the ground and the defined height and was" +
                   "discarded for gridding")
    return output_list_of_files

def create_gridded_dataset(Platform, Flightstr,Datestr, height_limit=10000, vertical_spacing=10):
    """
    Function to create gridded dataset from quality-controlled (Level-2) files
    if directory where NC files are stored is provided as a string, a list of file 
    paths for all NC files in the directory is created, otherwise a list of file 
    paths needed to be gridded can also be provided directly

    Arguments:
        Platform            (string):  Plane from which the dropsondes were released
        Flightstr           (string):  Flightnumber
        Datestr             (string):  Date of the flights
        height_limit        (integer): Height up to which the interpolation is performed 
                                       (default 10km)
        vertical_spacing    (integer): Vertical grid spacing (default 10m)                                                                                                                                                                                                                       
    Returns:
        dataset (xr.DataSet): Dataset with Level-3 structure
    """

    list_of_files  = sorted(glob.glob(archive + "Level_2/" + "*L2.nc"))

    list_of_files  = check_for_faults(list_of_files, height_limit, vertical_spacing)

    interp_list    = [None] * len(list_of_files)
    save_directory = archive + "Level_3/Interim_files/"
    
    if os.path.exists(archive + "Level_3") == False:                                        # check if directory already exists
            os.mkdir(archive + "Level_3/")                                                  # if not create it
    if os.path.exists(save_directory) == False:
            os.mkdir(save_directory)

    sonde_id = []

    for id_, file_path in enumerate(tqdm(list_of_files, desc=f"{Platform} {Flightstr}")):

        file_name = (str(file_path[-21:-5] + "L3.nc"))
        interp_list[id_] = interpolate_for_grid(file_path, 
                            height_limit=height_limit,
                            vertical_spacing=vertical_spacing)

        interp_list[id_].to_netcdf(save_directory + file_name, engine="netcdf4")
        sonde_id = np.append(sonde_id, file_path[-12:-6])

    concatenated_dataset = xr.concat(interp_list, dim="sounding")
    concatenated_dataset = concatenated_dataset.swap_dims({"sounding": "sonde_id"})
    concatenated_dataset = concatenated_dataset.assign_coords(sonde_id=sonde_id)

    return concatenated_dataset

#-------------------------------------------------------------------------------------------

for idx in np.arange(len(config["FLIGHT"]["Platforms"].split(", "))):                       # loop over all specified flights 
                                                                                            # in config file

    Platform  = config["FLIGHT"]["Platforms"].split(", ")[idx]                              # read out flight information
    Flightstr = config["FLIGHT"]["Flightstr"].split(", ")[idx]
    Datestr   = config["FLIGHT"]["Datestr"].split(", ")[idx]

    # define the directory where the data from the specified flight is stored in 
    # data_dir (defined in config), folder needs to have a folder called "Level_2" 
    # with *_L2.nc files from Level 2 processing
    archive = config["PATH"]["data_dir"] + config["FLIGHT"]["Campaignstr"] + \
              f"_{Platform}_Dropsondes_{Datestr}_{Flightstr}/"

    save_directory = archive + "/Level_3/"
    if os.path.exists(save_directory) == False:                                             # check if directory already exists,
        os.mkdir(save_directory)                                                            # if not create it
        print("Directory for saving Level 3 data created")

    gridded_dataset = create_gridded_dataset(Platform, Flightstr, Datestr,                  # gridding and interpolation of the data
                        height_limit=float(config["GRID SETTINGS"]["height_limit"]),
                        vertical_spacing=float(config["GRID SETTINGS"]["vertical_spacing"]))

    file_name = (f"{Datestr}_{Platform}_{Flightstr}_L3.nc")                                 # save gridded dataset to netcdf file
    gridded_dataset.to_netcdf(save_directory + file_name, mode="w", 
                              format="NETCDF4", engine="netcdf4")