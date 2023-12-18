import os
from configparser import ConfigParser
from importlib import reload

# read config file
config = ConfigParser(inline_comment_prefixes="#")
config.read("config.cfg")

if config["LEVELS"]["run_Level_2"] == "True":
    print("\nStarting Level-2 processing")
    if 'L2' in dir():
        reload(L2)
    else: 
        import Processing_scripts.Level_2 as L2
    print("Finished Level-2 processing")

if config["LEVELS"]["run_Level_3"] == "True":
    print("\nStarting Level-3 processing")
    if 'L3' in dir():
        reload(L3)
    else: 
        import Processing_scripts.Level_3 as L3
    print("Finished Level-3 processing")

if config["LEVELS"]["run_Level_4"] == "True":
    print("\nStarting Level-4 processing")
    if 'L4' in dir():
        reload(L4)
    else: 
        import Processing_scripts.Level_4 as L4
    print("Finished Level-4 processing")

if (config["LEVELS"]["Quicklooks_L3"]=="True") | (config["LEVELS"]["Quicklooks_L4"]=="True"):
    print("\nPlotting Quicklooks")
    if 'QL' in dir():
        reload(QL)
    else: 
        import Processing_scripts.Quicklooks as QL
    print("Finished Quicklooks\n")