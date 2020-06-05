"""
Convert Paul Chote's leodb JSON ephemerides to tlemcee format
"""

import os
import json
import numpy as np
import argparse as ap
from astropy.table import Table
from astropy.time import Time

def argParse():
    """
    Argument parser settings

    Returns
    -------
    args : array-like
        Array of command line arguments
    """
    parser = ap.ArgumentParser()

    parser.add_argument('json_path',
                        help='path to leodb JSON ephemeris file',
                        type=str)

    parser.add_argument('out_path',
                        help='path to output ephemeris file',
                        type=str)

    return parser.parse_args()

if __name__ == "__main__":

    args = argParse()

    if os.path.exists(args.json_path):
        with open(args.json_path) as f:
            json_data = json.load(f)
    else:
        print('error: input JSON file not found')
    
    # extract ephemeris info from JSON file
    mjd, ra, dec = [], [], []
    for eph in json_data['position_fits']:
        mjd.append(eph[0])
        ra.append(eph[1])
        dec.append(eph[2])
    
    # convert mjd to isot format
    utc = Time(mjd, format='mjd').isot
    
    ephem = Table()
    ephem['UTC'] = utc
    ephem['RA'] = ra
    ephem['DEC'] = dec
    
    # until uncertainty info included
    ephem['RAERR'] = (14. / 3600) / np.cos(np.deg2rad(ephem['DEC']))
    ephem['DECERR'] = 14. / 3600
    
    ephem = ephem[np.argsort(ephem['UTC'])]
    
    ephem.write(args.out_path, 
                format='csv', 
                overwrite=True)
