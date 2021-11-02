#!/bin/bash
#/* ================================================================================ *\
#!
#! *
#! * This file is part of pyMARS, a python-based package for IACT event analysis
#! * of the MAGIC telescopes. It is distributed to you in the hope that
#! * it can be a useful and timesaving tool in analysing Data of
#! * imaging Cerenkov telescopes. It is distributed WITHOUT ANY WARRANTY.
#! *
#! * Permission to use, copy, modify and distribute this software and its
#! * documentation for any purpose is hereby granted without fee,
#! * provided that the above copyright notice appear in all copies and
#! * that both that copyright notice and this permission notice appear
#! * in supporting documentation. It is provided "as is" without express
#! * or implied warranty.
#! *
#!
#!
#!   Author(s): Tjark Miener, 09/2021 <mailto:tmiener@ucm.es>
#!
#!   Copyright: pyMARS Software Development, 2021
#!
#!
#\* ================================================================================ */

#######################################################################################
#
#   pymars-merger tool for merging CTLearn prediction files
#
#   For help: $ pymars-merger -h
#
#######################################################################################

import argparse
import os
import glob
import numpy as np
import pandas as pd

def main():

    parser = argparse.ArgumentParser(
        description=("pymars-merger tool for merging CTLearn prediction files."))
    parser.add_argument('--particle', '-p',
                        help='input directory for the particle classification',
                        default=None)
    parser.add_argument('--energy', '-e',
                        help='input directory for the energy regression',
                        default=None)
    parser.add_argument('--arrdir', '-a',
                        help='input directory for the arrival direction regression',
                        default=None)
    parser.add_argument('--pattern', '-k',
                        help='pattern to mask unwanted files for the data input directory',
                        default=["*"],
                        nargs='+')
    parser.add_argument('--output_dir', '-o',
                        help='path where to save generated files. By default, the current directory is used.',
                        default="./")

    args = parser.parse_args()
    
    abs_file_dir = os.path.abspath(args.particle)
    print(abs_file_dir)
    filelists = []
    for pattern in args.pattern:
        print(pattern)
        files = glob.glob(os.path.join(abs_file_dir, pattern))
        print(files)
        if not files: continue
        filelists.extend(files)

    print(filelists)
    parameter_list = ['hillas_intensity', 'hillas_length', 'hillas_phi', 'hillas_psi', 'hillas_r', 'hillas_width',  'hillas_skewness', 'hillas_x', 'hillas_y', 'leakage_intensity_1', 'leakage_intensity_2']

    for file in filelists:
        filename = file.split("/")[-1]
        particle_data = pd.HDFStore(f"{file}", 'r')
        energy_data = pd.HDFStore(f"{args.energy}/{filename}", 'r')
        arrdir_data = pd.HDFStore(f"{args.arrdir}/{filename}", 'r')

        data = {}
        particle_prediction = particle_data['data']
        energy_prediction = energy_data['data']
        arrdir_prediction = arrdir_data['data']

        parameters = np.dstack(particle_prediction["parameters"].to_numpy())
        for i, parameter in enumerate(parameter_list):
            data["M1_"+parameter] = parameters[0][i]
            data["M2_"+parameter] = parameters[1][i]

        data['hadronness'] = 1.0 - np.array(particle_prediction['reco_gammaness']).astype(float)
        data['reco_alt'] = np.array(arrdir_prediction['reco_altitude']).astype(float)
        data['mc_alt'] = np.array(arrdir_prediction['mc_altitude']).astype(float)
        data['reco_az'] = np.array(arrdir_prediction['reco_azimuth']).astype(float)
        data['mc_az'] = np.array(arrdir_prediction['mc_azimuth']).astype(float)
        data['reco_energy'] = np.array(energy_prediction['reco_energy']).astype(float)
        data['pointing_alt'] = np.array(arrdir_prediction['pointing_alt']).astype(float)
        data['pointing_az'] = np.array(arrdir_prediction['pointing_az']).astype(float)

        if 'mc_particle' not in particle_prediction:
            data['mjd'] = np.array(particle_prediction['core_x']).astype(float)
            data['millisec'] = np.array(particle_prediction['core_y']).astype(float)
            data['nanosec'] = np.array(particle_prediction['mc_energy']).astype(float)
        else:
            data['mc_particle'] = np.array(particle_prediction['mc_particle']).astype(int)
            data['mc_energy'] = np.array(energy_prediction['mc_energy']).astype(float)

        if "data" in list(pd.HDFStore(f"{args.output_dir}/{filename}").keys()):
            pd.HDFStore(f"{args.output_dir}/{filename}").remove("data")
        pd.DataFrame(data=data).to_hdf(f"{args.output_dir}/{filename}", key="data", mode="a")

if __name__ == "__main__":
    main()
