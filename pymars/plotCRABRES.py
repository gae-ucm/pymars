
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
#   plotCRABRES - Python ctaplot-based program to plot angular resolution obtained
#                 with real observational data from MARS melibea and CTLearn files
#
#   For help: $ pymars-plotCRABRES -h
#
#######################################################################################

import argparse
import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy import units as u
import ctaplot
import uproot

def main():

    parser = argparse.ArgumentParser(
        description=("plotCRABRES - Python ctaplot-based program to plot angular resolution obtained with real observational data from MARS melibea and CTLearn files."))
    parser.add_argument('--input_data', '-i',
                        help='input directory',
                        default="./")
    parser.add_argument('--pattern', '-p',
                        help='pattern to mask unwanted files for the data input directory',
                        default="*")
    parser.add_argument('--input_mc', '-m',
                        help='input directory',
                        default=None)
    parser.add_argument('--pattern_mc', '-q',
                        help='pattern to mask unwanted files for the MCs input directory',
                        default="*")
    parser.add_argument('--output_dir', '-o',
                        help='path where to save generated files. By default, the current directory is used.',
                        default="./")
    parser.add_argument('--signalCut', '-s',
                        help='maximum theta2 values for the signal region',
                        default=0.025,
                        type=float)
    parser.add_argument('--hadronnessCut', '-k',
                        help='maximum hadronness values',
                        default=0.5,
                        type=float)
    parser.add_argument('--sizeCut', '-z',
                        help='minimum size values',
                        default=[50.0, 50.0],
                        nargs='+',
                        type=float)
    parser.add_argument('--leakage1Cut', '-l',
                        help='maximum leakage1 values',
                        default=[0.15, 0.15],
                        nargs='+',
                        type=float)
    parser.add_argument('--aleksic',
                        dest='aleksic',
                        help='plot MAGIC performance of Aleksic et al. (2015)',
                        default=False,
                        action='store_true')

    args = parser.parse_args()

    # Input handling
    abs_file_dir = os.path.abspath(args.input_data)
    filelists = [glob.glob(os.path.join(abs_file_dir, args.pattern))]
    if args.input_mc:
        abs_file_dir_mc = os.path.abspath(args.input_mc)
        filelists.append(glob.glob(os.path.join(abs_file_dir_mc, args.pattern_mc)))
        
    fig, ax = plt.subplots(1,figsize=(17,9))
    for f, files in enumerate(filelists):
        image_type = files[0].split("/")[-1].split("_")[1]
        label = f"Crab 2.93h ({image_type} images)" if f == 0 else f"MC ST0310 ({image_type} images)"
        for i, file in enumerate(files):
            print(file)
            data = pd.HDFStore(file, 'r')
            observations = data.keys() if f == 0 else ['gamma']
            for j, obs in enumerate(observations):
                reco_energy = np.array(data[obs]['reco_energy']) * u.TeV
                mc_alt = np.array(data[obs]['mc_altitude']) * u.rad
                reco_alt = np.array(data[obs]['reco_altitude']) * u.rad
                mc_az = np.array(data[obs]['mc_azimuth']) * u.rad
                reco_az = np.array(data[obs]['reco_azimuth']) * u.rad

                print(len(reco_energy))

                if args.signalCut or args.hadronnessCut or args.sizeCut or args.leakage1Cut:
                    parameters = np.dstack(data[obs]["parameters"].to_numpy())
                    hadronness = 1 - data[obs]["reco_gammaness"].to_numpy()
                    #TODO: Make generic for any given number of telescopes
                    size_m1 = parameters[0][0]
                    size_m2 = parameters[1][0]
                    leakage1_m1 = parameters[0][-2]
                    leakage1_m2 = parameters[1][-2]

                    if args.sizeCut:
                        #TODO: Make generic for any given number of telescopes
                        sizeMask = (size_m1 > args.sizeCut[0]) & (size_m2 > args.sizeCut[1])
                        reco_energy = reco_energy[sizeMask]
                        mc_alt = mc_alt[sizeMask]
                        reco_alt = reco_alt[sizeMask]
                        mc_az = mc_az[sizeMask]
                        reco_az = reco_az[sizeMask]
                        hadronness = hadronness[sizeMask]
                        leakage1_m1 = leakage1_m1[sizeMask]
                        leakage1_m2 = leakage1_m2[sizeMask]

                    if args.leakage1Cut:
                        #TODO: Make generic for any given number of telescopes
                        leakage1Mask = (leakage1_m1 < args.leakage1Cut[0]) & (leakage1_m2 < args.leakage1Cut[1])
                        reco_energy = reco_energy[leakage1Mask]
                        mc_alt = mc_alt[leakage1Mask]
                        reco_alt = reco_alt[leakage1Mask]
                        mc_az = mc_az[leakage1Mask]
                        reco_az = reco_az[leakage1Mask]
                        hadronness = hadronness[leakage1Mask]

                    if args.signalCut:
                        theta2 = (mc_alt-reco_alt).to(u.deg)**2 + (mc_az-reco_az).to(u.deg)**2
                        theta2Mask = (theta2.value < args.signalCut)
                        reco_energy = reco_energy[theta2Mask]
                        mc_alt = mc_alt[theta2Mask]
                        reco_alt = reco_alt[theta2Mask]
                        mc_az = mc_az[theta2Mask]
                        reco_az = reco_az[theta2Mask]
                        hadronness = hadronness[theta2Mask]

                    if args.hadronnessCut:
                        hadronnessMask = (hadronness < args.hadronnessCut)
                        reco_energy = reco_energy[hadronnessMask]
                        mc_alt = mc_alt[hadronnessMask]
                        reco_alt = reco_alt[hadronnessMask]
                        mc_az = mc_az[hadronnessMask]
                        reco_az = reco_az[hadronnessMask]

                print(len(reco_energy))

                if i == 0 and j == 0:
                    total_mc_alt = mc_alt
                    total_reco_alt = reco_alt
                    total_mc_az = mc_az
                    total_reco_az = reco_az
                    total_reco_energy = reco_energy
                else:
                    total_mc_alt = np.concatenate((total_mc_alt, mc_alt))
                    total_reco_alt = np.concatenate((total_reco_alt, reco_alt))
                    total_mc_az = np.concatenate((total_mc_az, mc_az))
                    total_reco_az = np.concatenate((total_reco_az, reco_az))
                    total_reco_energy = np.concatenate((total_reco_energy, reco_energy))

            print(len(total_reco_energy))

        ax = ctaplot.plot_angular_resolution_per_energy(total_reco_alt.to(u.rad).value, total_reco_az.to(u.rad).value, total_mc_alt.to(u.rad).value, total_mc_az.to(u.rad).value, total_reco_energy.to(u.TeV).value, bias_correction=False, label=f"{label}")

    if args.aleksic:
        # Aleksic et al. (2015) MAGIC performance (angular)
        aleksic_energy = np.array([80.27961169117373, 98.4965576498376, 149.39707865342058, 235.34800775765504, 379.26901907322497, 597.4695679970707, 948.3609962468593, 1493.9707865342073, 2389.4033799062745, 3792.6901907322535, 6020.121593466931, 9555.714328759786])*1e-3
        aleksic_resolution = np.array([0.16343127141077546, 0.15305784802242292, 0.12609730224229213, 0.10244695966988478, 0.08252238399252569, 0.06881286982248522, 0.0613156921519776, 0.05547600046714421, 0.05419115929616944, 0.05249279819370914, 0.05079443709124884, 0.054480691373403944])
        ax.plot(aleksic_energy, aleksic_resolution, label="Aleksic et al. (2015)")


    # Scale, labels and title
    ax.set_ylabel(r'$\theta [deg]$',fontsize=25)

    ax.set_ybound(0.03,0.12)
    ax.set_yticks([0.05,0.1])

    ax.set_xlabel("E_reco [TeV]",fontsize=25)
    ax.set_xscale('log')
    ax.set_xbound(4e-2,25)
    ax.set_title('Angular resolution',fontsize=30)

    ax.tick_params(labelsize=25)

    #major and minor grid lines
    plt.grid(b=True, which='major', color='gray', alpha=0.3, linestyle='dashdot', lw=1.5)
    plt.minorticks_on()
    plt.grid(b=True, which='minor', color='beige', alpha=0.5, ls='-', lw=1)
    ax.legend(loc='upper right',fontsize=25)

    plt.savefig(f"{args.output_dir}/angularresolution_crab.png")

if __name__ == "__main__":
    main()


