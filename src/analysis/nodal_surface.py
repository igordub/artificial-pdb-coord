# -*- coding: utf-8 -*-
import click
import logging
from pathlib import Path
from dotenv import find_dotenv, load_dotenv

from os import path
import glob
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import numpy as np
from pymol import cmd
from yaml import nodes
import src.utilities as utils
from src.analysis.modevectors import modevectors
import math

config = utils.read_config()
plt.style.use(config['viz']['default'])

eigenvector_filepaths = sorted(glob.glob("data/processed/eigenvector_???"))
no_modes = len(eigenvector_filepaths)
eigenvector_filepath_format = "data/processed/eigenvector_{:03}"
eigenvectors = {mode_num : np.loadtxt(eigenvector_filepath_format.format(mode_num)) for mode_num in range(7,no_modes+7)}

# Renormalize eigenvectors
# for key in eigenvectors.keys():
#     eigenvector = eigenvectors[key]
#     for row_idx in range(eigenvector.shape[0]):
#         eigenvector[row_idx] = eigenvector[row_idx] / np.linalg.norm(eigenvector[row_idx])



def calcualte_dotprod(eigenvector):
    no_beads = np.shape(eigenvector)[0]
    dotprod = np.zeros((no_beads, no_beads))

    for i in np.arange(no_beads):
        for j in np.arange(i+1):
            vector_1 = eigenvector[i]
            vector_2 = eigenvector[j]
            dotprod[i][j] = np.dot(vector_1, vector_2)

    dotprod = dotprod + dotprod.T

    return dotprod

def calcualte_crossprod_norm(eigenvector):
    no_beads = np.shape(eigenvector)[0]
    crossprod_norm = np.zeros((no_beads, no_beads))

    for i in np.arange(no_beads):
        for j in np.arange(i+1):
            vector_1 = eigenvector[i]
            vector_2 = eigenvector[j]
            crossprod_norm[i][j] = np.linalg.norm(np.cross(vector_1, vector_2))

    crossprod_norm = crossprod_norm + crossprod_norm.T

    return crossprod_norm

def plot_crosscor(array, vmin=0, center=None, vmax=1, axis=None, cmap=plt.cm.RdBu_r, annot=True):
#     array = array[1:, :-1]
#     mask = np.zeros_like(array)
#     mask[np.triu_indices_from(mask, k=1)] = True

    sns.heatmap(array, fmt='.2f', linewidths=.5, vmin=vmin, center=center, vmax=vmax, square=True, mask=None, annot=annot, cmap=cmap, cbar=False, ax=axis)
    
#     y_ticks = axis.get_yticks()
#     y_ticks = [int(tick +0.5) for tick in y_ticks]
#     axis.set_yticklabels(y_ticks)

    return None

def factorize(number):
    """ Factorizes number into two near integers.
    """
    x = math.floor(math.sqrt(number))
    y = math.ceil(number / x)
    return (x, y)

subplot_rows, subplot_columns = factorize(no_modes)

# with plt.rc_context({'font.size': 4, 'xtick.labelsize': 5, 'ytick.labelsize': 5, 'axes.titlesize': 7}):
_, axs = plt.subplots(subplot_rows, subplot_columns)
mode_num = 7
for i in range(subplot_rows):
    for j in range(subplot_columns):
        dotprod = calcualte_dotprod(eigenvectors[mode_num])
        plot_crosscor(dotprod, vmin=-1, center=0, vmax=1, cmap=plt.cm.RdBu_r, axis=axs[i, j], annot=False)

        axs[i, j].set_title("mode {}".format(mode_num))
        mode_num += 1

plt.tight_layout()
plt.savefig("scratch/dotprod.png")


# with plt.rc_context({'font.size': 4, 'xtick.labelsize': 5, 'ytick.labelsize': 5, 'axes.titlesize': 7}):
_, axs = plt.subplots(subplot_rows, subplot_columns)
mode_num = 7
for i in range(subplot_rows):
    for j in range(subplot_columns):
        corossprod_norm = calcualte_crossprod_norm(eigenvectors[mode_num])
        plot_crosscor(corossprod_norm, vmin=0, center=None, vmax=1, cmap=plt.cm.binary, axis=axs[i, j], annot=False)

        axs[i, j].set_title("mode {}".format(mode_num))
        mode_num += 1

plt.tight_layout()
plt.savefig("scratch/crossprod_norm.png")


# @click.command()
# @click.argument('input_dir', type=click.Path(exists=True))
# @click.argument('output_dir', type=click.Path())
# def main_commandline(input_dir, output_dir):
#     """ Runs data visualization scripts to turn processed data (from data/processed)
#         into plots (saved in scratch/).
#     """
#     logger = logging.getLogger(__name__)
#     logger.info('visualize results')

#     config = utils.read_config()
#     pdb_codes = config['pdb']['codeList']

#     plt.style.use(config['viz']['default'])
    
#     # Get filepaths
#     coords_filepaths = sorted(glob.glob(path.join(input_dir, "*.pdb")))


#     # Draw ENM
#     script_filepath = path.join(input_dir, "draw_enm.pml")
#     structure_filepath = path.join(input_dir, "CAonly.pdb")
#     cmd.reinitialize()
#     structure_name ="CAonly"
#     view = """set_view (\
#      0.254109710,    0.685759246,    0.682028472,\
#     -0.659577310,   -0.392885894,    0.640778780,\
#      0.707380295,   -0.612679243,    0.352475613,\
#      0.000000893,   -0.000001701,  -20.039798737,\
#      2.060348749,    0.924048603,    0.750656426,\
#     17.039802551,   23.039796829,  -20.000000000 )"""

#     draw_enm(structure_filepath, script_filepath, structure_name=structure_name, view=view)
    
#     cmd.set('sphere_scale', 0.5, structure_name)

#     cmd.viewport(width=1200, height=1200)
#     cmd.zoom(complete=1)

#     png_filepath = path.join(output_dir, "enm") + ".png"
#     pse_filepath = path.join(output_dir, "enm") + ".pse"

#     cmd.save(pse_filepath)
#     cmd.set('ray_opaque_background', 0)
#     cmd.png(png_filepath, width=1200, height=1200, ray=1)

#     # Draw eigenvectors
#     # cmd.reinitialize()
#     eigenvectors = {}
#     for mode_number in range(7,13):
#         first_structure = "CAonly"
#         mode_filepath = path.join(input_dir, "Mode_{:03}.pdb".format(mode_number))
#         last_structure = "struct_shift_{:03}".format(mode_number)
#         output_name= "mode_{:03}".format(mode_number)

#         cmd.load(mode_filepath, last_structure)
#         orange_colour = np.array([int(code) for code in config['colors']['orange'].split(',')])

#         eigenvector = modevectors(first_structure, last_structure, \
#             outname=output_name, cutoff=0, cut=0, factor=7.5, head=0.5, tail=0.2,\
#              head_length= 1, head_rgb = orange_colour/255, tail_rgb = orange_colour/255)     
#         cmd.delete(last_structure)

#         # Combine eigenvector components
#         eigenvectors[mode_number] = np.array(eigenvector).T
    
#     # Save eigenvectors
#     for mode_number, eigenvector in eigenvectors.items():
#         filename = "eigenvector_{:03}".format(mode_number)
#         np.savetxt(path.join("data/processed", filename), eigenvector, fmt='%.4e')

#     cmd.show_as('spheres', "CAonly")
#     cmd.show('sticks', "CAonly")
#     blue_colour = np.array([int(code) for code in config['colors']['blue'].split(',')])
#     cmd.color("0x{:02x}{:02x}{:02x}".format(*blue_colour), "rep spheres and CAonly")
#     cmd.set('grid_mode', 1)
#     cmd.set('grid_slot', -2, "CAonly")
#     cmd.set_view(\
#     (0.254109710,    0.685759246,    0.682028472,\
#     -0.659577310,   -0.392885894,    0.640778780,\
#      0.707380295,   -0.612679243,    0.352475613,\
#      0.000000000,    0.000000000,  -30.692567825,\
#      2.000000000,    1.154749870,    0.816750050,\
#    -79.449722290,  140.834869385,  -20.000000000))

#     cmd.save(path.join(output_dir, "modevectors.pse"))

# def main(input_dir, output_dir):
#     """ Runs data visualization scripts to turn processed data (from data/processed)
#         into plots (saved in scratch/).
#     """
#     logger = logging.getLogger(__name__)
#     logger.info('visualize results')

#     config = utils.read_config()
#     pdb_codes = config['pdb']['codeList']

#     mpl.rcParams.update(mpl.rcParamsDefault)  # VS Code plots not black
#     plt.style.use(config['viz']['default'])
    
#     # Get filepaths
#     coords_filepaths = sorted(glob.glob(path.join(input_dir, "Mode_???.pdb")))

#     # Draw ENM
#     script_filepath = path.join(input_dir, "draw_enm.pml")
#     structure_filepath = path.join(input_dir, "CAonly.pdb")
#     cmd.reinitialize()
#     structure_name ="CAonly"
#     pymol_view =    ((\
#                     0.254109710,    0.685759246,    0.682028472,\
#                     -0.659577310,   -0.392885894,    0.640778780,\
#                     0.707380295,   -0.612679243,    0.352475613,\
#                     0.000000893,   -0.000001701,  -20.039798737,\
#                     2.060348749,    0.924048603,    0.750656426,\
#                     17.039802551,   23.039796829,  -20.000000000 ))

#     draw_enm(structure_filepath, script_filepath, structure_name=structure_name, view=pymol_view)
    
#     cmd.set('sphere_scale', 0.5, structure_name)

#     cmd.viewport(width=1200, height=1200)
#     cmd.zoom(buffer=1,complete=1)


#     cmd.save(path.join(output_dir, "enm.pse"))
#     cmd.set('ray_opaque_background', 0)
#     cmd.png(path.join(output_dir, "enm.png"), width=1200, height=1200, ray=1)

#     # Draw eigenvectors
#     eigenvectors = {}
#     no_modes = len(coords_filepaths)
#     for mode_number in range(7, no_modes+1):
#         print("Mode number = {}".format(mode_number))
#         first_structure = "CAonly"
#         mode_filepath = path.join(input_dir, "Mode_{:03}.pdb".format(mode_number))
#         last_structure = "struct_shift_{:03}".format(mode_number)
#         output_name= "mode_{:03}".format(mode_number)

#         cmd.load(mode_filepath, last_structure)
#         orange_colour = np.array([int(code) for code in config['colors']['orange'].split(',')])

#         eigenvector = modevectors(first_structure, last_structure, \
#             outname=output_name, cutoff=0, cut=0, factor=5, head=0.5, tail=0.2,\
#              head_length= 1, head_rgb = orange_colour/255, tail_rgb = orange_colour/255)     
#         cmd.delete(last_structure)

#         # Combine eigenvector components
#         eigenvectors[mode_number] = np.array(eigenvector).T
    
#     # Save eigenvectors
#     for mode_number, eigenvector in eigenvectors.items():
#         filename = "eigenvector_{:03}".format(mode_number)
#         np.savetxt(path.join("data/processed", filename), eigenvector, fmt='%.4e')

#     cmd.show_as('spheres', "CAonly")
#     cmd.show('sticks', "CAonly")
#     blue_colour = np.array([int(code) for code in config['colors']['blue'].split(',')])
#     cmd.color("0x{:02x}{:02x}{:02x}".format(*blue_colour), "rep spheres and CAonly")
#     cmd.set('grid_mode', 1)
#     cmd.set('grid_slot', -2, "CAonly")

#     cmd.viewport(width=1600, height=800)
#     cmd.zoom(buffer=-5,complete=1)

#     cmd.save(path.join(output_dir, "modevectors.pse"))
#     cmd.ray()
#     cmd.png(path.join(output_dir, "modevectors.png"))