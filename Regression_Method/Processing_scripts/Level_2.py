import warnings

import datetime
import glob
import os

import numpy as np
import pandas as pd
import xarray as xr
from tqdm import trange
from configparser import ConfigParser

config = ConfigParser(inline_comment_prefixes="#")                                      # read config file
config.read("config.cfg")

varname_L1 = config["LEVEL 2"]["varname_L1"].split(", ")                                # load variable names to process
varname_L2 = config["LEVEL 2"]["varname_L2"].split(", ")                                # from level 1 to level 2
assert len(varname_L1) == len(varname_L2), "Lengths of variable list for Level 1 and Level 2 don't match."

alt_var    = config["LEVEL 2"]["altitude_varname"]                                      # get variable names from config file
lat_var    = config["LEVEL 2"]["latitude_varname"]
lon_var    = config["LEVEL 2"]["longitude_varname"]
time_var   = config["LEVEL 2"]["time_varname"]
p_var      = config["LEVEL 2"]["pressure_varname"]
rh_var     = config["LEVEL 2"]["relative_humidity_varname"]
ta_var     = config["LEVEL 2"]["air_temperature_varname"]
dp_var     = config["LEVEL 2"]["dew_point_varname"]

nc_meta = {                                                                             # definition of attributes for 
    time_var: {                                                                         # temperature, altitude, pressure,
        "long_name": config["LEVEL 2"]["time_longname"]},                               # latitude, longitude and relative
    alt_var:  {                                                                         # humidity with converted units
        "long_name": config["LEVEL 2"]["altitude_longname"],
        "units"    : "m"},
    lat_var:  {
        "long_name": config["LEVEL 2"]["latitude_longname"],
        "units"    : "degree north"},
    lon_var:  {
        "long_name": config["LEVEL 2"]["longitude_longname"],
        "units"    : "degree east"},
    p_var:    {
        "long_name": config["LEVEL 2"]["pressure_longname"],
        "units"    : "Pa"},
    ta_var:   {
        "long_name": config["LEVEL 2"]["air_temperature_longname"],
        "units"    : "K"},
    dp_var:   {
        "long_name": config["LEVEL 2"]["dew_point_longname"],
        "units"    : "K"},
    rh_var:   {
        "long_name": config["LEVEL 2"]["relative_humidity_longname"],
        "units"    : "Pa/Pa"}}

def get_all_sondes_list(Platform, Flightstr):
    """
    Extract lists for sondes in given directory

    Arguments: 
        Platform       (string): Name of the plane the dropsondes were released from

    Returns:
        sonde_ds       (list):   List with individual datasets of all sondes from PQC 
                                 files
        directory      (string): Directory with all Level_1 Files
        file_time_str  (list):   List with alle flight time strings
        sonde_paths    (list):   List with paths to the individual sonde file
    """
    
    directory = config["PATH"]["data_dir"] + config["FLIGHT"]["Campaignstr"] + \
                f"_{Platform}_Dropsondes_{Datestr}_{Flightstr}/Level_1/"                # define directory name where all sonde
                                                                                        # files are present

    sonde_paths = sorted(glob.glob(directory + "*QC.nc"))                               # get paths to the individual sonde
                                                                                        # files

    file_time_str = [None] * len(sonde_paths)                                           # lists to store extracted sonde time 
    file_time     = [None] * len(sonde_paths)                                           # from file name as string and as time
    sonde_ds      = [None] * len(sonde_paths)                                           # lists to store individual datasets 
                                                                                        # of all sondes from PQC filess

    for i in range(len(sonde_paths)):
        file_time_str[i] = sonde_paths[i][-20:-5]
        file_time[i]     = np.datetime64(pd.to_datetime(file_time_str[i], 
                                         format="%Y%m%d_%H%M%S"), "s")
        sonde_ds[i]      = xr.open_dataset(sonde_paths[i])

    return sonde_ds, directory, file_time_str, file_time, sonde_paths

def create_variable(ds, vname, data, **kwargs):
    """
    Insert the data into a variable in an :class:`xr.Dataset`
    
    Arguments:
        ds    (xr.Dataset): dataset to insert variable into
        vname (string)    : variable name
        data  (np.ndarray): variable array

    Returns:
        vname (string)
    """
    attrs = nc_meta[vname].copy()
    dims = [time_var]
    v = xr.Variable(dims, data, attrs=attrs)
    ds[vname] = v

    return vname

def get_global_attrs(Platform, file_time, sonde_ds):
    """
    Global attributes for output dataset

    Arguments: 
        Platform  (string):        Plane from which the dropsondes were released
        file_time (np.datetime64): release time
        sonde_ds  (xr.Dataset):    Level 1 processed dropsonde dataset

    Returns:
        nc_global_attrs (dict)
    """

    if hasattr(sonde_ds, "AspenVersion") == True:
            ASPENversion = sonde_ds.AspenVersion,
    elif hasattr(sonde_ds, "AvapsEditorVersion") == True:
            ASPENversion = sonde_ds.AvapsEditorVersion,

    nc_global_attrs = {
        "title": "Dropsondes processing Level-2",
        "platform_id":           Platform,
        "instrument_id":         config["INSTRUMENT"]["instrument_id"],
        "product_id":            "Level-2",
        "ASPEN_version":         ASPENversion,
        "ASPEN_processing_time": sonde_ds.ProcessingTime,
        "launch_date":           str(pd.to_datetime(file_time).date()),
        "author":                config["AUTHOR"]["name"],
        "author_email":          config["AUTHOR"]["email"],
        "author_institute":      config["AUTHOR"]["institute"],
        "creation_time":         str(datetime.datetime.utcnow()) + " UTC"
    }

    return nc_global_attrs
                                                                                        
#---------------------------------------------------------------------------------------

for idx in np.arange(len(config["FLIGHT"]["Platforms"].split(", "))):

    Platform  = config["FLIGHT"]["Platforms"].split(", ")[idx]
    Flightstr = config["FLIGHT"]["Flightstr"].split(", ")[idx]
    Datestr   = config["FLIGHT"]["Datestr"].split(", ")[idx]

    (sonde_ds, 
     directory, 
     file_time_str, 
     file_time, 
     sonde_paths) = get_all_sondes_list(Platform, Flightstr)

    save_directory = config["PATH"]["data_dir"] + config["FLIGHT"]["Campaignstr"] + \
                     f"_{Platform}_Dropsondes_{Datestr}_{Flightstr}/Level_2/"
    
    if os.path.exists(save_directory) == False:                                         # check if directory already exists,
        os.mkdir(save_directory)                                                        # if not create it
        print("Directory for saving Level 2 data created")

    for i in trange(len(sonde_ds), desc=f"{Platform} {Flightstr}"):

        ht_indices = ~np.isnan(sonde_ds[i][alt_var])                                    # retrieving non-NaN indices of 
                                                                                        # geopotential height (sonde_ds[i].gpsalt)
                                                                                        # only time values at these indices will 
                                                                                        # be used in Level-2 trajectory data

        obs = np.arange(1, ht_indices.sum() + 1, 1)                                     # creating the observations dimension
                                                                                        # of the NC file
        height              = sonde_ds[i][alt_var].values[ht_indices]                   # Variable array: geopotential height
        time                = sonde_ds[i][time_var].values[ht_indices]                  # Variable array: time
        variables           = {}
        variables[time_var] = time
        variables[alt_var]  = height

        if rh_var in varname_L1:                                                        # unit conversions of specified variables
            vname = varname_L2[varname_L1.index(rh_var)]
            variables[vname]  = np.float32(                                             # conversion of relative humidity
                                    sonde_ds[i][rh_var].values[ht_indices] / 100)       # from percent to fraction
        if p_var in varname_L1:
            vname = varname_L2[varname_L1.index(p_var)]
            variables[vname]  = np.float32(                                             # conversion of pressure
                                    sonde_ds[i][p_var].values[ht_indices] * 100)        # from hectopascal to pascal
        if ta_var in varname_L1:
            vname = varname_L2[varname_L1.index(ta_var)]                                                                                 
            variables[vname]  = np.float32(                                             # conversion of temperature
                                    sonde_ds[i][ta_var].values[ht_indices] +            # from celsius to kelvin
                                    float(config["CONSTANTS"]["tmelt"]))                        
        if dp_var in varname_L1:
            vname = varname_L2[varname_L1.index(dp_var)]
            variables[dp_var] = np.float32(                                             # conversion of dew point 
                                    sonde_ds[i][dp_var].values[ht_indices] +            # temperature from
                                    float(config["CONSTANTS"]["tmelt"]))                # celsius to kelvin

        for var1, var2 in zip(varname_L1, varname_L2):                                  # adding additional variables 
            if var2 not in variables.keys():                                            # without unit conversion
                variables[var2] = np.float32(sonde_ds[i][var1].values[ht_indices])      # add additional variables
            if var2 not in nc_meta.keys():                                              
                nc_meta[var2]   = sonde_ds[i][var1].attrs                               # and corresponding attributes

        to_save_ds = xr.Dataset(coords={time_var: obs})                                 # creating and populating dataset
        for var in nc_meta.keys():
            create_variable(to_save_ds, var, variables[var])

        sonde_id = file_time_str[i][-6:]                                                # adding a sonde_id var to the dataset
        attrs = {
            "descripion": "unique sonde ID",
            "long_name" : "sonde identifier",
        }
        sonde_id_var = xr.Variable([], sonde_id, attrs=attrs)
        to_save_ds["sonde_id"] = sonde_id_var

        file_name      =  file_time_str[i] + "_L2.nc"                                   # file name

        comp = dict(zlib=True, complevel=4, fletcher32=True,                            
                    _FillValue=np.finfo("float32").max)
        encoding = {var: comp for var in to_save_ds.data_vars if var != "sonde_id"}
        encoding["time"] = {"units": "seconds since 2020-01-01", "dtype": "float"}

        nc_global_attrs = get_global_attrs(Platform, file_time[i], sonde_ds[i])         # add global attributes
        for key in nc_global_attrs.keys():
            to_save_ds.attrs[key] = nc_global_attrs[key]

        to_save_ds.to_netcdf(save_directory + file_name, mode="w",                      # saving dataset to NetCDF file
                             format="NETCDF4", encoding=encoding)
