#!/usr/bin/env python3
import pkg_resources
from piranha.utils.log_colours import green,cyan
import sys
import os
from piranha.utils.config import *


def package_data_check(filename,directory,key,config):
    try:
        package_datafile = os.path.join(directory,filename)
        data = pkg_resources.resource_filename('piranha', package_datafile)
        config[key] = data
    except:
        sys.stderr.write(cyan(f'Error: Missing package data.')+f'\n\t- {filename}\nPlease install the latest piranha version with `piranha --update`.\n')
        sys.exit(-1)

def get_snakefile(thisdir,analysis_mode):

    snakefile = os.path.join(thisdir, 'scripts',f'piranha_{analysis_mode}.smk')
    if not os.path.exists(snakefile):
        sys.stderr.write(cyan(f'Error: cannot find Snakefile at {snakefile}\n Check installation or specify another analysis mode\n'))
        sys.exit(-1)
    return snakefile

def check_install(config):

    for resource in resources:
        package_data_check(resource[RESOURCE_KEY_FILENAME],resource[RESOURCE_KEY_DIRECTORY],resource[RESOURCE_KEY],config)

    if config[KEY_ANALYSIS_MODE] == "vp1":
        package_data_check(REFERENCE_SEQUENCES_FILE_VP1,"data",KEY_REFERENCE_SEQUENCES,config)
    elif "wg" in config[KEY_ANALYSIS_MODE]:
        package_data_check(REFERENCE_SEQUENCES_FILE_WG,"data",KEY_REFERENCE_SEQUENCES,config)

# config={}
# check_install()
