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
def main(output_dir):
    """ Creates artificial PDB files (saved in pdb/raw).
    """
    logger = logging.getLogger(__name__)
    logger.info('creating artificial PDB files')

    config = utils.read_config()

    tetrahed_coords = make_tetrahedron(bond_length=10)
    # print(tetrahed_coords)
    write_pdb_file(tetrahed_coords, output_filepath=path.join(output_dir, "tetrahed.pdb"))

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

    main()
