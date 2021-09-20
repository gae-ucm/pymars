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
#   podie - Python odie-like program to generate the theta2 plot automatically from
#           MARS melibea (under development) and CTLearn files.
#
#   For help: $ pymars-podie -h
#
#######################################################################################

import argparse
import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy import units as u

def main():

    parser = argparse.ArgumentParser(
        description=("podie - Python odie-like program to generate the theta2 plot automatically from MARS melibea (under development) and CTLearn files."))
    parser.add_argument('--input_dir', '-i',
                        help='input directory',
                        default="./")
    parser.add_argument('--pattern', '-p',
                        help='pattern to mask unwanted files',
                        default="*.h5")
    parser.add_argument('--output_dir', '-o',
                        help='path where to save generated files. By default, the current directory is used.',
                        default="./")
    parser.add_argument('--name', '-n',
                        help='name or ID of the anaylsis',
                        default=None)
    parser.add_argument('--eRange', '-e',
                        help='target energy range; valid option are (FR, LE, HE)',
                        default=None)
    parser.add_argument('--signalCut', '-s',
                        help='maximum theta2 values for the signal region (overwrites default value of eRange argument)',
                        default=None,
                        type=float)
    parser.add_argument('--hadronnessCut', '-k',
                        help='maximum hadronness values (overwrites default value of eRange argument)',
                        default=None,
                        type=float)
    parser.add_argument('--sizeCut', '-z',
                        help='minimum size values (automatically set if eRange provided)',
                        default=[300.0, 300.0],
                        nargs='+',
                        type=float)
    parser.add_argument('--energyCut',
                        help='minimum reco energy values (automatically set if eRange provided)',
                        default=0.0,
                        type=float)
    parser.add_argument('--nWobbleOff', '-w',
                        help='number of off regions',
                        default=None,
                        type=int)
    parser.add_argument('--offRegions',
                        help='list of off regions with rotation angle in degrees (automatically set if nWobbleOff provided; default are 3 off regions with "[180, 90, 270]" degrees)',
                        default=[180, 90, 270],
                        nargs='+',
                        type=int)
    parser.add_argument('--rangeTh2', '-r',
                        help='range of the theta2 plot; default "0.4"',
                        default=0.4,
                        type=float)
    parser.add_argument('--nBinsSignal', '-b',
                        help='number of bins in signal region; default "2"',
                        default=2,
                        type=int)

    args = parser.parse_args()

    # Input handling
    abs_file_dir = os.path.abspath(args.input_dir)
    files = glob.glob(os.path.join(abs_file_dir, args.pattern))
    filename_type = files[0].split('.')[-1]

    # MARS standard analysis
    # For full range (FR) analysis
    # Sensitivity ~ 0.7% Crab
    # For low energy (LE) analysis
    # Sensitivity ~ 1.2% Crab
    # For high energy (HE) analysis
    # Sensitivity ~ 1.% Crab
    eRanges = {
        "FR": {
            "signalCut": 0.009,
            "hadronnessCut": 0.16,
            "sizeCut": [300.0, 300.0],
            "energyCut": 0.0
        },
        "LE": {
            "signalCut": 0.02,
            "hadronnessCut": 0.28,
            "sizeCut": [60.0, 60.0],
            "energyCut": 0.0
        },
        "HE": {
            "signalCut": 0.007,
            "hadronnessCut": 0.1,
            "sizeCut": [400.0, 400.0],
            "energyCut": 1000.0
        },
    }
    
    # Parse the quality cuts for this analysis
    eRange = {}
    if args.eRange:
        eRange = eRanges[args.eRange]
    else:
        eRange["sizeCut"] = args.sizeCut
        eRange["energyCut"] = args.energyCut
    if args.hadronnessCut: eRange["hadronnessCut"] = args.hadronnessCut
    if args.signalCut: eRange["signalCut"] = args.signalCut
    
    # Set up the histogram bins
    hist_bin = int(args.nBinsSignal*args.rangeTh2/eRange["signalCut"])
    
    # Set up the off regions
    nWobbleOff2rotations = {1: [180], 3:[180, 90, 270]}
    offRegions = args.offRegions
    if args.nWobbleOff: offRegions = nWobbleOff2rotations[args.nWobbleOff]
    
    if filename_type == "h5":

        off_managment = {}
        for i, file in enumerate(files):
            data = pd.HDFStore(file, 'r')
            observations = data.keys()
            for j, obs in enumerate(observations):
                off_managment[obs] = pd.DataFrame({})
                
                parameters = np.dstack(data[obs]["parameters"].to_numpy())
                
                #TODO: Make generic for any given number of telescopes
                size_m1 = parameters[0][0]
                size_m2 = parameters[1][0]
                sizeMask = (size_m1 > eRange["sizeCut"][0]) & (size_m2 > eRange["sizeCut"][1])
                
                off_managment[obs] = pd.DataFrame({})

                hadronness = 1 - data[obs]["reco_gammaness"].to_numpy()

                hadronness_mask = (hadronness[sizeMask] < eRange["hadronnessCut"])
                reco_energy = data[obs]["reco_energy"].to_numpy()
                energyMask = (reco_energy[sizeMask][hadronness_mask] > eRange["energyCut"])
                # ON regions
                theta2_on = ((data[obs]["mc_altitude"]-data[obs]["reco_altitude"]).to_numpy() * u.rad).to(u.deg)**2 + ((data[obs]["mc_azimuth"]-data[obs]["reco_azimuth"]).to_numpy() * u.rad).to(u.deg)**2
                theta2_on_hadronness = theta2_on[sizeMask][hadronness_mask][energyMask]
                if i == 0 and j == 0:
                    total_on_hadronness = theta2_on_hadronness.value
                else:
                    total_on_hadronness = np.concatenate((total_on_hadronness, theta2_on_hadronness.value))

                # OFF regions
                delta_alt = (data[obs]["mc_altitude"]-data[obs]["pointing_alt"]).to_numpy()
                delta_az = (data[obs]["mc_azimuth"]-data[obs]["pointing_az"]).to_numpy()
                delta = np.vstack((delta_alt, delta_az))
                for rotation in offRegions:
                    rotation = float(rotation) * u.deg
                    rotation_matrix = np.matrix([[np.cos(rotation.to(u.rad)), -np.sin(rotation.to(u.rad))],
                                             [np.sin(rotation.to(u.rad)), np.cos(rotation.to(u.rad))]], dtype=float)
                    off = np.squeeze(np.asarray(np.dot(rotation_matrix, delta)))
                    off_managment[obs][f"off_alt_{rotation.value}"] = off[:][0]
                    off_managment[obs][f"off_az_{rotation.value}"] = off[:][1]
                    off_managment[obs][f"off_alt_{rotation.value}_pointing"] = off_managment[obs][f"off_alt_{rotation.value}"] + data[obs]["pointing_alt"]
                    off_managment[obs][f"off_az_{rotation.value}_pointing"] = off_managment[obs][f"off_az_{rotation.value}"] + data[obs]["pointing_az"]

                
        # Create the theta2 plot
        fig, ax = plt.subplots(figsize =(10, 7))

        # Plot the histogram for the On region
        ax.hist(total_on_hadronness, bins=hist_bin, range=(0.0, args.rangeTh2), color= "black", alpha=1.0, histtype='step', label='On')

        # Plot the histogram for the Off regions
        n_off = 0
        for rotation in offRegions:
            rotation = float(rotation)
            for i, file in enumerate(files):
                data = pd.HDFStore(file, 'r')
                observations = data.keys()
                for j, obs in enumerate(observations):
                   
                    parameters =  np.dstack(data[obs]["parameters"].to_numpy())
                    #TODO: Make generic for any given number of telescopes
                    size_m1 = parameters[0][0]
                    size_m2 = parameters[1][0]
                    sizeMask = (size_m1 > eRange["sizeCut"][0]) & (size_m2 > eRange["sizeCut"][1])
                    hadronness = 1 - data[obs]["reco_gammaness"].to_numpy()
                    hadronness_mask = (hadronness[sizeMask] < eRange["hadronnessCut"])
                    
                    reco_energy = data[obs]["reco_energy"].to_numpy()
                    energyMask = (reco_energy[sizeMask][hadronness_mask] > eRange["energyCut"])
                    theta2_off = ((off_managment[obs][f"off_alt_{rotation}_pointing"]-data[obs]["reco_altitude"]).to_numpy() * u.rad).to(u.deg)**2 + ((off_managment[obs][f"off_az_{rotation}_pointing"]-data[obs]["reco_azimuth"]).to_numpy() * u.rad).to(u.deg)**2
                    theta2_off_hadronness = theta2_off[sizeMask][hadronness_mask][energyMask]

                    if i == 0 and j == 0:
                        total_off_hadronness = theta2_off_hadronness.value
                    else:
                        total_off_hadronness = np.concatenate((total_off_hadronness, theta2_off_hadronness.value))
                    
            ax.hist(total_off_hadronness, bins=hist_bin, range=(0.0, args.rangeTh2), alpha=0.5, histtype='step', label=f'Off_{int(rotation)}')
            n_off += len(total_off_hadronness[total_off_hadronness < eRange["signalCut"]])
            
        n_on = float(len(total_on_hadronness[total_on_hadronness < eRange["signalCut"]]))
        n_off = float(n_off)
        n_signal = n_on + n_off
        alpha = 1.0/float(len(offRegions))
        log_on = np.log( (n_on*(alpha + 1)) / (n_signal*alpha))
        log_off = np.log( (n_off*(alpha + 1)) / n_signal)
        sig = np.sqrt((n_on*log_on + n_off*log_off) * 2)

        gamma_rate = (n_on - n_off/float(len(offRegions)))/(2.93*60)
        bkg_rate = (n_off/float(len(offRegions)))/(2.93*60)

        ax.text(0.8, 0.9, f'Time = 2.93 h', fontsize=12, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
        ax.text(0.8, 0.85, f'N_on = {int(n_on)}; N_off = {int(n_off/float(len(offRegions)))}', fontsize=12, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
        ax.text(0.8, 0.8, f'Significance (Li&Ma) = {sig:.1f}$ \sigma $', fontsize=12, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
        ax.text(0.8, 0.75, f'Gamma Rate = {gamma_rate:.2f} / min', fontsize=12, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
        ax.text(0.8, 0.7, f'Bkg Rate = {bkg_rate:.3f} / min', fontsize=12, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

        ax.set_xbound(0.0, args.rangeTh2)
        ax.axvline(x=eRange["signalCut"], linestyle="dashed", color="red")
        ax.legend(loc='upper center',fontsize=17)

        ax.text(0.5, 0.5, 'Preliminary', fontsize=58, alpha=0.5, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

        ax.set_title(f"{args.name}",fontsize=25)
        ax.set_xlabel(r'$\theta^{2} [deg^{2}]$', fontsize=25)
        ax.set_ylabel(r'$N_{events}$', fontsize=25)
        ax.tick_params(labelsize=20)

        # Save plot
        plt.savefig(f"{args.output_dir}/podie_{args.name}.png")
        # Show plot
        plt.show()

if __name__ == "__main__":
    main()


