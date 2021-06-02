#!/usr/bin/env python
""" This is the master script for recreating the results.

    It imports each of the key module fucntions 
    and runs them one by one.

    Run the whole thing from the root directory 
    to replicate all results:
    
    $ python -m src.main
"""

import src.utilities as utils

config = utils.read_config()

import src.structure.make_pdb as mkpdb
import src.simulation.simulate_enm as simenm
import src.analysis.visualize as visual

mkpdb.main_import(config['pdb']['rawFilePath'])
simenm.main_import(config['pdb']['rawFilePath'], config['data']['rawFilePath'])
visual.main_import(config['data']['rawFilePath'], config['data']['outPathScratch'])
