import warnings

import datetime
import os
import globe
from tqdm import tqdm

import xarray as xr
import numpy as np
from configparser import ConfigParser

import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

config = ConfigParser(inline_comment_prefixes="#")											# read config file
config.read("config.cfg")

alt_var    = config["LEVEL 3"]["altitude_varname"]

def plot_quicklooks_L3():

	directory = archive + "/Level_3/"

	if os.path.exists(directory) == False:                                        			# check if directory already exists
		print("Level 3 directory not found")

	file_name = (f"{Datestr}_{Platform}_{Flightstr}_L3.nc")
	dataset = xr.open_dataset(directory + file_name)

	n_subplots = len(list(dataset.variables.keys())) - len(									# get number of variables 
				 config["QUICKLOOKS"]["skip_variables_L3"].split(", "))						# (number of subplots)

	col_subplots = int(config["QUICKLOOKS"]["number_of_columns"])							# get number of columns
	row_subplots = int(np.ceil(n_subplots/col_subplots))									# calculate number of rows

	fig, ax = plt.subplots(nrows=row_subplots, ncols=col_subplots, 							# create figure with subplots
						   figsize=(col_subplots*2, row_subplots*3))

	col_count, row_count = 0, 0  															# create counters for subplots
	
	for par in dataset.variables:															# loop over all variables to plot
	
		if par not in config["QUICKLOOKS"]["skip_variables_L3"].split(", "):				# leave out variables defined in 
																							# Default section of config	

			if "sonde_id" in dataset[par].dims:												# for all variables with data for 
																							# each sonde, plot mean as black
				av = dataset[par].mean(dim="sonde_id")										# line and individual profiles 
				ax[row_count, col_count].plot(av.data, av[alt_var], color="black", 			# as light grey lines
											label="mean")
				for idx in dataset.sonde_id.data:
					sonde = dataset[par].sel(sonde_id=idx)
					ax[row_count, col_count].plot(sonde.data, sonde[alt_var], 
											color="grey", alpha=0.4)

			else:						
				ax[row_count, col_count].plot(dataset[par].data, dataset[alt_var], 			# for all other variables plot
											color="black")									# profiles in black

			if col_count != 0: ax[row_count, col_count].set_yticklabels([])					# add labels and axes description
			else: ax[row_count, col_count].set_ylabel("height [m]")							# only add ylabels for first plot row
			ax[row_count, col_count].set_xlabel(f"{par} [{dataset[par].attrs['units']}]")

			if col_count != 3: col_count = col_count + 1 									# add to the row and column counter
			else: col_count, row_count = 0, row_count + 1

	fig.suptitle(f"{Datestr} {Platform} {Flightstr}", fontsize=16)							# add title to plot
 
	while col_count <= 3: 																	# delete unused axes
		fig.delaxes(ax[row_count][col_count])
		col_count = col_count+1

	plt.tight_layout()

	if os.path.exists("Quicklooks/") == False: os.mkdir("Quicklooks/")						# save plot in Quicklooks folder
	plt.savefig(f"Quicklooks/L3_quicklook_{Datestr}_{Platform}_{Flightstr}.png", dpi=300)
	plt.close()

	return

def plot_quicklooks_L4():

	directory = archive + "/Level_4/"

	if os.path.exists(directory) == False:                                        			# check if directory already exists
		print("Level 4 directory not found")

	if config["REGRESSION"]["coordinate_system"] == "carthesian":							# define coordinate names depending
		coor1 = "x"																			# on coordinate system used in
		coor2 = "y"																			# Level 4 processing
	elif config["REGRESSION"]["coordinate_system"] == "spherical":
		coor1 = "theta"
		coor2 = "phi"
	else: print("Please enter coordinate system option, either spherical or carthesian")

	for pattern_id in config["PATTERN IDS"][Flightstr].split(", "):							# loop over all flight patterns

		file_name = (f"{Datestr}_{Platform}_{Flightstr}_{pattern_id}_L4.nc")
		dataset = xr.open_dataset(directory + file_name)									# get dataset from file

		n_subplots = len(config["QUICKLOOKS"]["plot_singlevar_L4"].split(", ")) + \
					 len(config["QUICKLOOKS"]["plot_gradients_L4"].split(", "))				# get number of variables 
					 																		# (number of subplots)

		col_subplots = int(config["QUICKLOOKS"]["number_of_columns"])						# get number of subplot columns
		row_subplots = int(np.ceil(n_subplots/col_subplots))								# calculate number of subplot rows

		fig, ax = plt.subplots(nrows=row_subplots, ncols=col_subplots, 						# create figure with subplots
							   figsize=(col_subplots*2, row_subplots*3))

		col_count, row_count = 0, 0 														# define counters for subplots

		for par in config["QUICKLOOKS"]["plot_gradients_L4"].split(", "):					# loop over all gradient variables
																							# defined in Default section of config

			ax[row_count, col_count].plot(dataset[f"d{par}_d{coor1}"].data, 				# plot profile for first coordinate
												   dataset[alt_var], color="black",			# in black
										           linestyle="dashed", label=coor1)
			ax[row_count, col_count].fill_betweenx(dataset[alt_var].data,			 		# plot standard error of profile as
												   dataset[f"d{par}_d{coor1}"].data - 		# grey area
												   dataset[f"d{par}_d{coor1}_stderr"].data,
												   dataset[f"d{par}_d{coor1}"].data +
												   dataset[f"d{par}_d{coor1}_stderr"].data,
												   color="grey", alpha=0.3,)		

			ax[row_count, col_count].plot(dataset[f"d{par}_d{coor2}"].data, 				# plot profile for second coordinate 
												   dataset[alt_var], color="black",
												   label=coor2)
			ax[row_count, col_count].fill_betweenx(dataset[alt_var].data, 					# plot standard error of profile
												   dataset[f"d{par}_d{coor2}"].data -
											       dataset[f"d{par}_d{coor2}_stderr"].data, 
												   dataset[f"d{par}_d{coor2}"].data + 
												   dataset[f"d{par}_d{coor2}_stderr"].data,
												   color="grey", alpha=0.3)
			ax[row_count, col_count].legend()

			if col_count != 0: ax[row_count, col_count].set_yticklabels([])					# add labels and axes description
			else: ax[row_count, col_count].set_ylabel("height [m]")							# only add ylabel and ticks for
			ax[row_count, col_count].set_xlabel(f"{par} [{dataset[par].attrs['units']}]")	# first plot in row

			if col_count != 3: col_count = col_count + 1 									# add to row and column counter
			else: col_count, row_count = 0, row_count + 1
			
		for par in config["QUICKLOOKS"]["plot_singlevar_L4"].split(", "):					# loop over all single variables
																							# defined in Default section of config

			ax[row_count, col_count].plot(dataset[par], dataset[alt_var], color="black")	# plot profile in black
			ax[row_count, col_count].fill_betweenx(dataset[alt_var],						# add standard error as grey area
												   dataset[par] - dataset[f"{par}_stderr"],
												   dataset[par] + dataset[f"{par}_stderr"],
												   color="grey", alpha=0.3)

			if col_count != 0: ax[row_count, col_count].set_yticklabels([])					# add labels and axes description
			else: ax[row_count, col_count].set_ylabel("height [m]")							# only add ylabel and ticks for 
			ax[row_count, col_count].set_xlabel(f"{par} [{dataset[par].attrs['units']}]")	# first plot in row

			if col_count != 3: col_count = col_count + 1 									# add to row and column counter
			else: col_count, row_count = 0, row_count + 1

		fig.suptitle(f"{Datestr} {Platform} {Flightstr} {pattern_id}", fontsize=16)			# add title to plot
	 
		while col_count <= 3: 																# delete all empty subplots
			fig.delaxes(ax[row_count][col_count])
			col_count = col_count+1

		plt.tight_layout()

		if os.path.exists("Quicklooks/") == False: os.mkdir("Quicklooks/")					# save plot to Quicklooks folder
		plt.savefig(f"Quicklooks/L4_quicklook_{Datestr}_{Platform}_{Flightstr}_{pattern_id}.png")
		plt.close()

	return

for idx in np.arange(len(config["FLIGHT"]["Platforms"].split(", "))):                       # loop over all specified flights 
                                                                                            # in config file
    
    Platform  = config["FLIGHT"]["Platforms"].split(", ")[idx]								# read out flight information
    Flightstr = config["FLIGHT"]["Flightstr"].split(", ")[idx]
    Datestr   = config["FLIGHT"]["Datestr"].split(", ")[idx]

    archive = config["PATH"]["data_dir"] + config["FLIGHT"]["Campaignstr"] + \
              f"_{Platform}_Dropsondes_{Datestr}_{Flightstr}/"

    if config["LEVELS"]["Quicklooks_L3"] == "True":
    	plot_quicklooks_L3()

    if config["LEVELS"]["Quicklooks_L4"] == "True":
   		plot_quicklooks_L4()
