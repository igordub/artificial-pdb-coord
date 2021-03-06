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

# utils.clean()
# mkpdb.main(config['pdb']['rawFilePath'])
# simenm.main(config['pdb']['rawFilePath'], config['data']['rawFilePath'])
visual.main(config['data']['rawFilePath'], config['data']['outPathScratch'])
