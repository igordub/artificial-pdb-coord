# -*- coding: utf-8 -*-
import click
import logging
from pathlib import Path
from dotenv import find_dotenv, load_dotenv

import subprocess
from os import path
import glob
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import numpy as np
from pymol import cmd
import src.utilities as utils
from src.analysis.modevectors import modevectors
import src.analysis.viz_1point as viz_1point


@click.command()
@click.argument('input_dir', type=click.Path(exists=True))
@click.argument('output_dir', type=click.Path())
def main(input_dir, output_dir):
    """ Runs data visualization scripts to turn processed data (from data/processed)
        into plots (saved in scratch/).
    """
    logger = logging.getLogger(__name__)
    logger.info('visualize results')

    config = utils.read_config()
    pdb_codes = config['pdb']['codeList']

    mpl.rcParams.update(mpl.rcParamsDefault)  # VS Code plots not black
    plt.style.use(config['viz'])
    
    # Get filepaths
    coords_filepaths = sorted(glob.glob(path.join(input_dir, "*.pdb")))


    # Draw ENM
    script_filepath = path.join(input_dir, "draw_enm.pml")
    structure_filepath = path.join(input_dir, "CAonly.pdb")
    cmd.reinitialize()
    structure_name ="tetrahed"
    view = """set_view (\
     0.254109710,    0.685759246,    0.682028472,\
    -0.659577310,   -0.392885894,    0.640778780,\
     0.707380295,   -0.612679243,    0.352475613,\
     0.000000893,   -0.000001701,  -20.039798737,\
     2.060348749,    0.924048603,    0.750656426,\
    17.039802551,   23.039796829,  -20.000000000 )"""
    draw_enm(structure_filepath, script_filepath, structure_name=structure_name, \
                output_dir=output_dir, view=view)
    
    # Draw eigenvectors
    # cmd.reinitialize()
    for mode_number in range(7,13):
        first_structure = path.join(input_dir, "CAonly.pdb")
        last_structure = path.join(input_dir, "Mode_{:03}.pdb".format(mode_number))
        output_name = "mode_{}".format(mode_number)
        draw_eigenvecotrs(first_structure, last_structure, output_name)
    cmd.save(path.join(output_dir, "modevectors.pse"))

def load_data(data_filepath):
    """ Load interim data into dataframe.
    """
    # Load data
    interim_data = pd.read_csv(data_filepath, sep=',', header=0)

    return interim_data

def draw_enm(structure_filepath, script_filepath, structure_name="enm", output_dir='.', view=None):
    """ Draws elastic network mode_numberl of a structure and save image.
    """

    cmd.delete('all')
    cmd.load(structure_filepath, "CAonly")
    cmd.run(script_filepath)

    # Set name
    if structure_name == "CAonly":
        pass
    else:
        cmd.set_name("CAonly", structure_name)

    # Set view
    if view == None:
        cmd.orient()
    else:
        cmd.set_view(view)
    
    cmd.set('sphere_scale', 0.5, structure_name)
    cmd.viewport(width=1200, height=1200)
    cmd.zoom(complete=1)

    png_filepath = path.join(output_dir, structure_name) + ".png"
    pse_filepath = path.join(output_dir, structure_name) + ".pse"

    cmd.save(pse_filepath)
    cmd.set('ray_opaque_background', 0)
    cmd.png(png_filepath, width=1200, height=1200, ray=1)

    return (pse_filepath, png_filepath)

def draw_eigenvecotrs(first_structure, last_structure, output_name):
    """ Draws eginvectors based on two PDB structure coordinate difference.
    """
    cmd.delete("first_obj_frame")
    cmd.delete("last_obj_frame")
    cmd.delete(output_name)
    cmd.load(first_structure, object="first_obj_frame")
    cmd.load(last_structure, object="last_obj_frame")

    modevectors("first_obj_frame", "last_obj_frame", outname=output_name,cutoff=0, cut=0, factor=5, head=0.5, tail=0.2)

    cmd.delete("first_obj_frame")
    cmd.delete("last_obj_frame")

if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    # not used in this stub but often useful for finding various files
    project_dir = Path(__file__).resolve().parents[2]

    # find .env automagically by walking up directories until it's found, then
    # load up the .env entries as environment variables
    load_dotenv(find_dotenv())

    main()
