# -*- coding: utf-8 -*-
import click
import logging
from pathlib import Path
from dotenv import find_dotenv, load_dotenv

from os import path
import numpy as np
import src.utilities as utils
from urllib.request import urlretrieve

@click.command()
@click.argument('output_dir', type=click.Path())
def main_commandline(output_dir):
    """ Creates artificial PDB files (saved in pdb/raw).
    """
    logger = logging.getLogger(__name__)
    logger.info('creating artificial PDB files')

    config = utils.read_config()

    cube_coords = make_cube(bond_length=4)
    # print(tetrahed_coords)
    write_pdb_file(cube_coords, output_filepath=path.join(output_dir, "struct.0.pdb"))

    return None

def main(output_dir):
    """ Creates artificial PDB files (saved in pdb/raw).
    """
    logger = logging.getLogger(__name__)
    logger.info('creating artificial PDB files')

    config = utils.read_config()

    coords = make_tetrahedron(bond_length=4)
    # coords = make_cube(bond_length=4)

    write_pdb_file(coords, output_filepath=path.join(output_dir, "struct.0.pdb"))

    return None

def make_tetrahedron(bond_length=4):
    """ Creates coordinates of an tetrahedron, equilateral triangular pyramid.
        Returns coordinate array.
        
        Input arguments:
        bond_length - measured in angstroms
    """
    vertices = np.array([   (0,0,0),
                            (1,0,0),
                            (0.5, np.sqrt(3)*0.5,0),
                            (0.5,np.sqrt(3)/6,np.sqrt(2/3))])
    coords = bond_length * vertices

    return coords

def make_cube(bond_length=4):
    """ Creates coordinates of a cube.
        Returns coordinate array.
        
        Input arguments:
        bond_length - measured in angstroms
    """
    vertices = np.array([   (0,0,0),
                            (1,0,0),
                            (0,1,0),
                            (0,0,1),
                            (1,1,0),
                            (1,0,1),
                            (0,1,1),
                            (1,1,1)])
    coords = bond_length * vertices

    return coords

def make_plane(no_x, no_y, bond_length=4):
    """ Creates coordinates of a plane.
        Returns coordinate array.
        
        Input arguments:
        bond_length - measured in angstroms
    """
    x_values = np.arange(no_x)
    y_values = np.arange(no_y)
    z_values = [0] 

    # x_coords, y_coords = np.meshgrid(x_values, y_values, sparse=True)
    # z_coords = [z_value]

    x_coords, y_coords, z_coords = (x_values, y_values, z_values)

    coords = []
    for x_coord in x_coords:
        for y_coord in y_coords:
            for z_coord in z_coords:
                coord = (x_coord, y_coord, z_coord)
                coords.append(coord)

    coords = np.array(coords) * bond_length

    return coords

def make_cuboid(no_x, no_y, no_z, bond_length=4):
    """ Creates coordinates of a plane.
        Returns coordinate array.
        
        Input arguments:
        bond_length - measured in angstroms
    """
    x_values = np.arange(no_x)
    y_values = np.arange(no_y)
    z_values = np.arange(no_z) 

    x_coords, y_coords, z_coords = (x_values, y_values, z_values)

    coords = []
    for x_coord in x_coords:
        for y_coord in y_coords:
            for z_coord in z_coords:
                coord = (x_coord, y_coord, z_coord)
                coords.append(coord)

    coords = np.array(coords) * bond_length

    return coords

def write_pdb_file(coords, output_filepath="structure.pdb"):
    """ Writes a PDB file from a coordinate array.
    """
    record_name = "ATOM"
    alt_loc = " "
    atom_name = "CA"
    residue_name = "ALA"
    chain_id = "A"
    insertion = " "
    occupancy = 1.00
    bfactor = 0.00
    segment_id = " "
    element_symbol = "C"
    charge = " "

    record_line = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n"

    with open(output_filepath, 'w') as file:
        for idx, coord in enumerate(coords):
            atom_number = idx + 1
            residue_number = idx + 1
            x_coord = coord[0]
            y_coord = coord[1]
            z_coord = coord[2]

            record = (record_name, atom_number, atom_name, alt_loc, residue_name, chain_id, \
                residue_number, insertion, x_coord, y_coord, z_coord, occupancy, bfactor, \
                segment_id, element_symbol, charge)
            
            file.write(record_line.format(*record))
    
    return None

        

def download_pdb(pdb_codes, output_filepath, biounit = True, compressed = False):
    """ Downloads raw PDB files form a list of PDB IDs.
        Authored by Chris Swain (http://www.macinchem.org)
        Modified by Igors Dubanevics (https://github.com/igordub)
        Copyright CC-BY
    """
    
    for pdb_code in pdb_codes:
        # Add .pdb extension and remove ':1' suffix in entities
        filename = "{:4s}.pdb".format(pdb_code[:4])
        
        # Add '1' if biounit
        if biounit:
            filename = "{}1".format(filename)
        # Add .gz extenison if compressed
        elif compressed:
            filename = "{}.gz".format(filename)
        
        url = path.join("https://files.rcsb.org/download/", filename.lower())
        destination_file = path.join(output_filepath, filename)
        # Download file
        urlretrieve(url, destination_file)

    return None


if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)

    # not used in this stub but often useful for finding various files
    project_dir = Path(__file__).resolve().parents[2]

    # find .env automagically by walking up directories until it's found, then
    # load up the .env entries as environment variables
    load_dotenv(find_dotenv())

    main_commandline()
