# python script to write all possible combinations of dropsondes to a .txt file
# output file contains a counter number followed by 4 spaces and the 
# list of sonde ids divided each by " ,"
# for each combination "run_processing.py" is called
#
# created 03.04.2023 by Fiona Paulus IGMK (fpaulus@uni-koeln.de)

from itertools import chain, combinations
import datetime
import os
import fileinput
import sys
from tqdm import tqdm
from importlib import reload
from csv import writer
from configparser import ConfigParser
import xarray as xr
import numpy as np

# Function definition ----------------------------------------------

def powerset(list_name, min_length):
    """
    Function to give out powerset fo list of all possible 
    combinations

    Arguments:
        list_name  (list):    name of the list to process
        min_length (integer): minimum number of entries per
                              combination
    
    Returns chain of possible combinations
    """
    s = list(list_name)                                             # make sure input is a list
    return chain.from_iterable(combinations(s, r) for r in range(   # create chain of possible 
                               min_length, len(s)+1))               # combinations

# Sonde ids and minimum length -------------------------------------

min_num_sondes = 5                                                  # minumum numboer of sondes 
                                                                    # for combination list

RF10_01 	= ["134740", "135441", "135946", "140443", "141152"]    # arrays with all sonde ids 
RF11_01 	= ["095337", "095819", "100245", "100704", "101127",    # for each flight with 
               "101559", "102036", "102530", "104454", "104822"]    # dropsonde pattern during
RF11_02 	= ["110301", "110700", "111058", "111511", "111938",    # HALO-(AC)3 campaign
               "112405", "112819", "113212", "114240", "114632"]
RF11_03 	= ["120144", "120537", "120926", "121319", "121736", 
               "122206", "122631", "123036", "124254", "124629"]
RF08_01 	= ["121910", "123313", "124741", "125824", "132412", 
               "133336", "140211", "141210", "142505", "143807", 
               "131122", "150552", "151423", "152724", "134817"] 
RF18_01		= ["093216", "093546", "093904", "094217", "094526", 
               "094841", "095157", "095517", "095842", "100214"]
RF15_01		= ["065242", "065754", "070314", "070904", "071549", 
               "072440", "073054", "073822", "074521", "075108", 
               "061904", "062708", "075441", "083047", "090409", 
               "091001", "100242"]

listofflights = [RF08_01, RF18_01] #[RF10_01, RF11_01, RF11_02, RF11_03, RF08_01,       # list of flights to process
                 #RF18_01]       
flightinfo    = [#"HALO-AC3 20220329 HALO RF10 C01",
                 #"HALO-AC3 20220330 HALO RF11 C01",                 # list of information about
                 #"HALO-AC3 20220330 HALO RF11 C02",                 # flights as headers for 
                 #"HALO-AC3 20220330 HALO RF11 C03",                 # output file
                 "HALO-AC3 20220330 P5   RF08 C01",
                 "HALO-AC3 20220412 HALO RF18 C01"]

os.system("cp config.cfg config_backup.cfg")                        # create backup of original
                                                                    # config file

# Specifications for outfile ---------------------------------------

author    = "Fiona Paulus"                                          # set author of output file
institute = "IGMK, Cologne"                                         # research institute
email     = "fpaulus@uni-koeln.de"                                  # contact email

"""with open('combinations.csv', 'w', newline='') as f:

    writer_object = writer(f)
    writer_object.writerow([f"All possible combinations of sondes "+ # write header with information
                            f"with minimum length of "             + # of file contents
                            f"{min_num_sondes}"])
    writer_object.writerow([f"created by {author}, {institute} "   + # write contact details to file
                            f"({email})"])
    writer_object.writerow([datetime.datetime.now().strftime(        # add creation date
                            "%d %B %Y %I:%M%p")])                 
    writer_object.writerow([])                                       # add empty line
    
    writer_object.writerow(["Flight", "ID", "number of sondes",      # add header for output
                            "sonde ids", "mean divergence [1/s]",
                            "mean divergence stderr [1/s]",
                            "mean vorticity [1/s]",
                            "mean vorticity stderr [1/s]"])
    f.close()"""

# Function calls ---------------------------------------------------

for idx, listofsondes in enumerate(listofflights):                  # loop over all flights in list     
    
    with open('combinations.csv', 'a', newline='') as f:            
        writer_object = writer(f)
        writer_object.writerow([f"{flightinfo[idx]}"])              # add flight information
        f.close()
    
    counter = 1

    for x in powerset(listofsondes, min_num_sondes):                # create chain of combinations
        
        str = x[0]                                                  # define output string with
                                                                    # first sonde id
        for idstr in x[1:]:                                         # loop over all chain entries
            str = str + ", " + idstr                                # and add sonde id to output

        config_orig = open('config_backup.cfg','r')                 # open backup config file
        config      = open('config.cfg','w')                        # open config file to alter

        for line in config_orig:                                    # loop over all lines in 
                                                                    # backup config file
        
            if line.startswith("Platforms"):                        # add platform to altered
                config.write(f"Platforms 	= "                   + # config file for running
                             f"{flightinfo[idx][18:22]} \n")        # combination
            elif line.startswith("Datestr"):                        # add datestring
                config.write(f"Datestr 	= "                       +
                             f"{flightinfo[idx][9:17]} \n")
            elif line.startswith("Flightstr"):                      # add flight string
                config.write(f"Flightstr 	= "                   +
                             f"{flightinfo[idx][23:27]} \n")
            elif line.startswith("[PATTERN IDS]"):                  # add pattern id and comment
                config.write(line)                                  
                config.write(f"{flightinfo[idx][23:27]} = "       +        
                             f"{flightinfo[idx][-2:]} \n")
            elif line.startswith("[PATTERN SONDE IDS]"):            # add sonde ids for combination
                config.write(line)
                config.write(f"{flightinfo[idx][23:27]}_"         +
                             f"{flightinfo[idx][-2:]} = {str} \n")
            else:
                config.write(line)                                  # all other lines are copied                                                        
                                                                    # directly

        config.close()                                              # close config file and
        config_orig.close()                                         # backup file

        print(f"\nprocessing {flightinfo[idx]} combination "      +
              f"{counter} ({len(x)} sondes)")
        
        configf = ConfigParser(inline_comment_prefixes="#")          # read config file
        configf.read("config.cfg")

        if flightinfo[idx][18:22] == "HALO":
            archive = configf["PATH"]["data_dir"]                 +\
                      f"{flightinfo[idx][:8]}_"                   +\
                      f"{flightinfo[idx][18:22]}_Dropsondes_"     +\
                      f"{flightinfo[idx][9:17]}_"                 +\
                      f"{flightinfo[idx][23:27]}/Level_4/"          # define archive path
            filename= f"{flightinfo[idx][9:17]}_"                 +\
                      f"{flightinfo[idx][18:22]}_"                +\
                      f"{flightinfo[idx][23:27]}_"                +\
                      f"{flightinfo[idx][-2:]}_L4.nc"               # define filename for HALO
        else:
            archive = configf["PATH"]["data_dir"]                 +\
                      f"{flightinfo[idx][:8]}_"                   +\
                      f"{flightinfo[idx][18:20]}_Dropsondes_"     +\
                      f"{flightinfo[idx][9:17]}_"                 +\
                      f"{flightinfo[idx][23:27]}/Level_4/"          # define archive path
            filename= f"{flightinfo[idx][9:17]}_"                 +\
                      f"{flightinfo[idx][18:20]}_"                +\
                      f"{flightinfo[idx][23:27]}_"                +\
                      f"{flightinfo[idx][-2:]}_L4.nc"               # define filename for P5
                
        if os.path.exists(archive + filename):                      # check if archive folder8
            os.system(f"rm {archive+filename}")                     # delete data file to overwrite

        if 'run_processing' in dir():                               # if it has already been imported
            reload(run_processing)
        else:
            import run_processing                                   # run processing script
       
        ds = xr.open_dataset(archive + filename)                    # load dataset to calculate mean
                                                                    # divergence and standard error

        nds = len(ds.sonde_id)                                      # number of sondes in dataset
    
        mask = (ds.div < 10**30) & (~np.isnan(ds.div_stderr.data))  # mask for valid data
        mean_div    = ds.div[mask].mean(dim="gpsalt").data          # calculate mean divergence
        mean_div_se = ds.div_stderr[mask].mean(dim="gpsalt").data   # calculate mean standard error
        mean_vor_se = ds.vor_stderr[mask].mean(dim="gpsalt").data   # calculate mean standard error
        mean_vor    = ds.vor[mask].mean(dim="gpsalt").data          # calculate mean vorticity

        del ds                                                      # delete dataset to overwrite
        
        with open('combinations.csv', 'a', newline='') as f:
            writer_object = writer(f)
            writer_object.writerow([flightinfo[idx][-8:], counter,  # write results to file
                                    nds, str, mean_div, mean_div_se, 
                                    mean_vor, mean_vor_se])
            f.close()
                                           
        counter = counter + 1

os.system("cp config_backup.cfg config.cfg")                        # create backup of original
                                                                    # config file
os.system("rm config_backup.cfg")                                   # delete backup of config file