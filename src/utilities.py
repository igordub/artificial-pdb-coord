#!/usr/bin/env python
"""
This script provides useful funcs to all other scripts
"""
import yaml
import os
import glob

def clean(flag="all"):
    """ Cleans project directories.
    """
    subdirs = ["raw", "interim", "processed"]
    data_filepaths = []
    pdb_filepaths = []
    for subdir in subdirs:
        data = glob.glob("data/{}/*".format(subdir))
        data_filepaths.extend(data)
        pdb = glob.glob("pdb/{}/*".format(subdir))
        pdb_filepaths.extend(pdb)

    scratch = glob.glob("scratch/*")

    if flag == "all":
        for filepath in data_filepaths:
            os.remove(filepath)
        for filepath in pdb_filepaths:
            os.remove(filepath)
        for filepath in scratch:
            os.remove(filepath)

    return None

def read_config():
    # Read in config file
    with open("config.yaml") as yaml_file:
        # YAML loads a list of dictionaries
        config_list = yaml.full_load(yaml_file)
        # Convert list into dict
        config_dict = {key: value for dict in config_list for key, value in dict.items()}
    return config_dict

def despine(ax, spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward', 10))  # outward by 10 points
        else:
            spine.set_color('none')  # don't draw spine

    # turn off ticks where there is no spine
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        # no yaxis ticks
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        # no xaxis ticks
        ax.xaxis.set_ticks([])