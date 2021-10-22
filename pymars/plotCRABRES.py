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
                        default=["*"],
                        nargs='+')
    parser.add_argument('--input_mc', '-m',
                        help='input directory',
                        default=None)
    parser.add_argument('--pattern_mc', '-q',
                        help='pattern to mask unwanted files for the MCs input directory',
                        default=["*"],
                        nargs='+')
    parser.add_argument('--output_dir', '-o',
                        help='path where to save generated files. By default, the current directory is used.',
                        default="./")
    parser.add_argument('--theta2Cut', '-t',
                        help='maximum theta2 values',
                        default=0.3,
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
    parser.add_argument('--nWobbleOff', '-w',
                        help='number of off regions',
                        default=None,
                        type=int)
    parser.add_argument('--offRegions',
                        help='list of off regions with rotation angle in degrees (automatically set if nWobbleOff provided; default are 3 off regions with "[180, 90, 270]" degrees)',
                        default=[180, 90, 270],
                        nargs='+',
                        type=int)
    parser.add_argument('--nBins', '-b',
                        help='number of bins; default "20"',
                        default=20,
                        type=int)
    parser.add_argument('--aleksic',
                        dest='aleksic',
                        help='plot MAGIC performance of Aleksic et al. (2015)',
                        default=False,
                        action='store_true')

    args = parser.parse_args()

    # Input handling
    abs_file_dir = os.path.abspath(args.input_data)
    filelists = []
    labels = []
    for pattern in args.pattern:
        files = glob.glob(os.path.join(abs_file_dir, pattern))
        if not files: continue
        filelists.append(files)
        if files[0].split('.')[-1] == "h5":
            labels.append(f"Crab 2.93h ({files[0].split('/')[-1].split('_')[1]} images)")
        elif files[0].split('.')[-1] == "root":
            labels.append(f"MARS - Crab 2.93h")

    if args.input_mc:
        abs_file_dir_mc = os.path.abspath(args.input_mc)
        for pattern in args.pattern_mc:
            files = glob.glob(os.path.join(abs_file_dir, pattern))
            if not files: continue
            filelists.append(files)
            if files[0].split('.')[-1] == "h5":
                labels.append(f"MC ST0310 ({files[0].split('/')[-1].split('_')[1]} images)")
            elif files[0].split('.')[-1] == "root":
                labels.append(f"MARS - MC ST0310")
    print(labels)

    # Set up the off regions
    nWobbleOff2rotations = {1: [180], 3:[180, 90, 270]}
    offRegions = args.offRegions
    if args.nWobbleOff: offRegions = nWobbleOff2rotations[args.nWobbleOff]

    fig, ax = plt.subplots(1,figsize=(17,9))
    for f, files in enumerate(filelists):
        if files[0].split('.')[-1] == "h5":
            for i, file in enumerate(files):
                print(file)
                data = pd.HDFStore(file, 'r')
                observations = data.keys() if  "Crab" in labels[f] else ['gamma']
                for j, obs in enumerate(observations):
                    reco_energy = np.array(data[obs]['reco_energy']) * u.TeV
                    mc_alt = np.array(data[obs]['mc_altitude']) * u.rad
                    reco_alt = np.array(data[obs]['reco_altitude']) * u.rad
                    mc_az = np.array(data[obs]['mc_azimuth']) * u.rad
                    reco_az = np.array(data[obs]['reco_azimuth']) * u.rad
                    pointing_alt = np.array(data[obs]['pointing_alt']) * u.rad
                    pointing_az = np.array(data[obs]['pointing_az']) * u.rad
                    print(len(reco_energy))

                    if args.hadronnessCut or args.sizeCut or args.leakage1Cut:
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
                            pointing_alt = pointing_alt[sizeMask]
                            pointing_az = pointing_az[sizeMask]
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
                            pointing_alt = pointing_alt[leakage1Mask]
                            pointing_az = pointing_az[leakage1Mask]
                            hadronness = hadronness[leakage1Mask]

                        if args.hadronnessCut:
                            hadronnessMask = (hadronness < args.hadronnessCut)
                            reco_energy = reco_energy[hadronnessMask]
                            mc_alt = mc_alt[hadronnessMask]
                            reco_alt = reco_alt[hadronnessMask]
                            mc_az = mc_az[hadronnessMask]
                            reco_az = reco_az[hadronnessMask]
                            pointing_alt = pointing_alt[hadronnessMask]
                            pointing_az = pointing_az[hadronnessMask]

                    # ON region
                    theta2_on = (mc_alt-reco_alt).to(u.deg)**2 + (mc_az-reco_az).to(u.deg)**2
                    if  "Crab" in labels[f]:
                        reco_energy_on = reco_energy
                        mc_alt_on = mc_alt
                        reco_alt_on = reco_alt
                        mc_az_on = mc_az
                        reco_az_on = reco_az
                        # OFF regions
                        #off_counts = np.zeros(args.nBins)
                        delta = np.vstack((mc_alt-pointing_alt, mc_az-pointing_az))
                        for r, rotation in enumerate(offRegions):
                            rotation = float(rotation) * u.deg
                            rotation_matrix = np.matrix([[np.cos(rotation.to(u.rad)), -np.sin(rotation.to(u.rad))],
                                                        [np.sin(rotation.to(u.rad)), np.cos(rotation.to(u.rad))]], dtype=float)
                            off = np.squeeze(np.asarray(np.dot(rotation_matrix, delta)))
                            print(off)
                            off_alt = off[:][0] * u.rad + pointing_alt
                            off_az = off[:][1] * u.rad + pointing_az
                            print(off_alt)
                            print(off_az)
                            print(pointing_alt)
                            print(pointing_az)
                            
                            theta2_off = (off_alt-reco_alt).to(u.deg)**2 + (off_az-reco_az).to(u.deg)**2
                            off_counts_rot, off_bins = np.histogram(theta2_off, bins=args.nBins, range=(0.0, args.theta2Cut))
                            if i == 0 and j == 0 and r == 0:
                                off_counts = off_counts_rot
                            else:
                                off_counts += off_counts_rot

                        print(off_counts)
                        print(off_bins)


                    else:
                        theta2Mask = (theta2_on.value < args.theta2Cut)
                        reco_energy_on = reco_energy[theta2Mask]
                        mc_alt_on = mc_alt[theta2Mask]
                        reco_alt_on = reco_alt[theta2Mask]
                        mc_az_on = mc_az[theta2Mask]
                        reco_az_on = reco_az[theta2Mask]
                    print(len(reco_energy_on))

                    if i == 0 and j == 0:
                        total_theta2_on = theta2_on
                        total_mc_alt = mc_alt_on
                        total_reco_alt = reco_alt_on
                        total_mc_az = mc_az_on
                        total_reco_az = reco_az_on
                        total_reco_energy = reco_energy_on
                    else:
                        total_theta2_on = np.concatenate((total_theta2_on, theta2_on))
                        total_mc_alt = np.concatenate((total_mc_alt, mc_alt_on))
                        total_reco_alt = np.concatenate((total_reco_alt, reco_alt_on))
                        total_mc_az = np.concatenate((total_mc_az, mc_az_on))
                        total_reco_az = np.concatenate((total_reco_az, reco_az_on))
                        total_reco_energy = np.concatenate((total_reco_energy, reco_energy_on))

                print(len(total_reco_energy))

        elif files[0].split('.')[-1] == "root":
            for i, file in enumerate(files):
                print(file)
                melibea_file = uproot.open(file)
                evts = melibea_file["Events"]
                reco_alt = np.asarray(evts["MStereoParDisp.fDirectionY"].array()) +  90.0 - np.asarray(evts["MPointingPos_1.fZd"].array())
                reco_az = np.asarray(evts["MStereoParDisp.fDirectionX"].array()) +  np.asarray(evts["MPointingPos_1.fAz"].array())
                marsDefaultMask = ~np.isnan(reco_alt)
                reco_alt *= u.deg
                reco_az *= u.deg
                mc_alt =  (np.asarray(evts["MSrcPosCam_1.fY"].array()) * 0.00337 +  90.0 - np.asarray(evts["MPointingPos_1.fZd"].array())) * u.deg
                mc_az =  (np.asarray(evts["MSrcPosCam_1.fX"].array()) * 0.00337 +  np.asarray(evts["MPointingPos_1.fAz"].array())) * u.deg
                pointing_alt = (90.0 - np.asarray(evts["MPointingPos_1.fZd"].array())) * u.deg
                pointing_az = np.asarray(evts["MPointingPos_1.fAz"].array()) * u.deg
                reco_energy = np.asarray(evts["MEnergyEst.fEnergy"].array()) * u.GeV


                print(len(reco_energy))

                if args.hadronnessCut or args.sizeCut or args.leakage1Cut:
                    reco_energy = reco_energy[marsDefaultMask]
                    mc_alt = mc_alt[marsDefaultMask]
                    reco_alt = reco_alt[marsDefaultMask]
                    mc_az = mc_az[marsDefaultMask]
                    reco_az = reco_az[marsDefaultMask]
                    pointing_alt = pointing_alt[marsDefaultMask]
                    pointing_az = pointing_az[marsDefaultMask]
                    hadronness = np.asarray(evts["MHadronness.fHadronness"].array())[marsDefaultMask]
                    size_m1 = np.asarray(evts["MHillas_1.fSize"].array())[marsDefaultMask]
                    size_m2 = np.asarray(evts["MHillas_2.fSize"].array())[marsDefaultMask]
                    leakage1_m1 = np.asarray(evts["MNewImagePar_1.fLeakage1"].array())[marsDefaultMask]
                    leakage1_m2 = np.asarray(evts["MNewImagePar_2.fLeakage1"].array())[marsDefaultMask]

                    if args.sizeCut:
                        #TODO: Make generic for any given number of telescopes
                        sizeMask = (size_m1 > args.sizeCut[0]) & (size_m2 > args.sizeCut[1])
                        reco_energy = reco_energy[sizeMask]
                        mc_alt = mc_alt[sizeMask]
                        reco_alt = reco_alt[sizeMask]
                        mc_az = mc_az[sizeMask]
                        reco_az = reco_az[sizeMask]
                        pointing_alt = pointing_alt[sizeMask]
                        pointing_az = pointing_az[sizeMask]
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
                        pointing_alt = pointing_alt[leakage1Mask]
                        pointing_az = pointing_az[leakage1Mask]
                        hadronness = hadronness[leakage1Mask]

                    if args.hadronnessCut:
                        hadronnessMask = (hadronness < args.hadronnessCut)
                        reco_energy = reco_energy[hadronnessMask]
                        mc_alt = mc_alt[hadronnessMask]
                        reco_alt = reco_alt[hadronnessMask]
                        mc_az = mc_az[hadronnessMask]
                        reco_az = reco_az[hadronnessMask]
                        pointing_alt = pointing_alt[hadronnessMask]
                        pointing_az = pointing_az[hadronnessMask]


                theta2_on = (mc_alt-reco_alt).to(u.deg)**2 + (mc_az-reco_az).to(u.deg)**2
                if  "Crab" in labels[f]:
                    reco_energy_on = reco_energy
                    mc_alt_on = mc_alt
                    reco_alt_on = reco_alt
                    mc_az_on = mc_az
                    reco_az_on = reco_az
                    # OFF regions
                    #off_counts = np.zeros(args.nBins)
                    delta = np.vstack(((mc_alt-pointing_alt).to(u.rad), (mc_az-pointing_az).to(u.rad)))
                    for r, rotation in enumerate(offRegions):
                        rotation = float(rotation) * u.deg
                        rotation_matrix = np.matrix([[np.cos(rotation.to(u.rad)), -np.sin(rotation.to(u.rad))],
                                                    [np.sin(rotation.to(u.rad)), np.cos(rotation.to(u.rad))]], dtype=float)
                        off = np.squeeze(np.asarray(np.dot(rotation_matrix, delta)))
                        print(off)
                        off_alt = off[:][0] * u.rad + pointing_alt
                        off_az = off[:][1] * u.rad + pointing_az
                        print(off_alt)
                        print(off_az)
                        print(pointing_alt)
                        print(pointing_az)

                        theta2_off = (off_alt-reco_alt).to(u.deg)**2 + (off_az-reco_az).to(u.deg)**2
                        off_counts_rot, off_bins = np.histogram(theta2_off, bins=args.nBins, range=(0.0, args.theta2Cut))
                        if i == 0 and r == 0:
                            off_counts = off_counts_rot
                        else:
                            off_counts += off_counts_rot

                    print(off_counts)
                    print(off_bins)

                else:
                    theta2Mask = (theta2_on.value < args.theta2Cut)
                    reco_energy_on = reco_energy[theta2Mask]
                    mc_alt_on = mc_alt[theta2Mask]
                    reco_alt_on = reco_alt[theta2Mask]
                    mc_az_on = mc_az[theta2Mask]
                    reco_az_on = reco_az[theta2Mask]
                print(len(reco_energy_on))

                if i == 0:
                    total_theta2_on = theta2_on
                    total_mc_alt = mc_alt_on
                    total_reco_alt = reco_alt_on
                    total_mc_az = mc_az_on
                    total_reco_az = reco_az_on
                    total_reco_energy = reco_energy_on
                else:
                    total_theta2_on = np.concatenate((total_theta2_on, theta2_on))
                    total_mc_alt = np.concatenate((total_mc_alt, mc_alt_on))
                    total_reco_alt = np.concatenate((total_reco_alt, reco_alt_on))
                    total_mc_az = np.concatenate((total_mc_az, mc_az_on))
                    total_reco_az = np.concatenate((total_reco_az, reco_az_on))
                    total_reco_energy = np.concatenate((total_reco_energy, reco_energy_on))


        if  "Crab" in labels[f]:

            on_counts, on_bins = np.histogram(total_theta2_on, bins=args.nBins, range=(0.0, args.theta2Cut))
            off_counts //= len(offRegions)
            print(off_counts)
            for b, bin in enumerate(off_bins):
                if b == 0:
                    previous_bin = bin
                    continue

                print(b)
                print(bin)
                print(previous_bin)
                off_counts_in_bin = int(off_counts[b-1])
                print(off_counts[b-1])
                theta2Mask = (total_theta2_on > previous_bin) & (total_theta2_on < bin)
                previous_bin = bin


                reco_energy_ex = total_reco_energy[theta2Mask]
                mc_alt_ex = total_mc_alt[theta2Mask]
                reco_alt_ex = total_reco_alt[theta2Mask]
                mc_az_ex = total_mc_az[theta2Mask]
                reco_az_ex = total_reco_az[theta2Mask]

                print(on_counts[b-1])
                print("alala")
                print(len(reco_energy_ex))
                if len(reco_energy_ex) > off_counts_in_bin:
                    if b == 1:
                        print("here1")
                        total_mc_alt_ex = mc_alt_ex[off_counts_in_bin:]
                        total_reco_alt_ex = reco_alt_ex[off_counts_in_bin:]
                        total_mc_az_ex = mc_az_ex[off_counts_in_bin:]
                        total_reco_az_ex = reco_az_ex[off_counts_in_bin:]
                        total_reco_energy_ex = reco_energy_ex[off_counts_in_bin:]
                    else:
                        print("here2")
                        total_mc_alt_ex = np.concatenate((total_mc_alt_ex, mc_alt_ex[off_counts_in_bin:]))
                        total_reco_alt_ex = np.concatenate((total_reco_alt_ex, reco_alt_ex[off_counts_in_bin:]))
                        total_mc_az_ex = np.concatenate((total_mc_az_ex, mc_az_ex[off_counts_in_bin:]))
                        total_reco_az_ex = np.concatenate((total_reco_az_ex, reco_az_ex[off_counts_in_bin:]))
                        total_reco_energy_ex = np.concatenate((total_reco_energy_ex, reco_energy_ex[off_counts_in_bin:]))

            ax = ctaplot.plot_angular_resolution_per_energy(total_reco_alt_ex.to(u.rad).value, total_reco_az_ex.to(u.rad).value, total_mc_alt_ex.to(u.rad).value, total_mc_az_ex.to(u.rad).value, total_reco_energy_ex.to(u.TeV).value, bias_correction=False, label=f"{labels[f]}")
        else:
            ax = ctaplot.plot_angular_resolution_per_energy(total_reco_alt.to(u.rad).value, total_reco_az.to(u.rad).value, total_mc_alt.to(u.rad).value, total_mc_az.to(u.rad).value, total_reco_energy.to(u.TeV).value, bias_correction=False, label=f"{labels[f]}")

    if args.aleksic:
        # Aleksic et al. (2015) MAGIC performance (angular)
        aleksic_energy = np.array([80.27961169117373, 98.4965576498376, 149.39707865342058, 235.34800775765504, 379.26901907322497, 597.4695679970707, 948.3609962468593, 1493.9707865342073, 2389.4033799062745, 3792.6901907322535, 6020.121593466931, 9555.714328759786])*1e-3
        aleksic_resolution = np.array([0.16343127141077546, 0.15305784802242292, 0.12609730224229213, 0.10244695966988478, 0.08252238399252569, 0.06881286982248522, 0.0613156921519776, 0.05547600046714421, 0.05419115929616944, 0.05249279819370914, 0.05079443709124884, 0.054480691373403944])
        ax.plot(aleksic_energy, aleksic_resolution, label="Aleksic et al. (2015)")

    # Scale, labels and title
    ax.set_ylabel(r'$\theta_{68} [deg]$',fontsize=25)

    ax.set_ybound(0.03,0.12)
    ax.set_yticks([0.05,0.1])

    ax.set_xlabel("Reconstructed energy [TeV]",fontsize=25)
    ax.set_xscale('log')
    ax.set_xbound(0.1,10)
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


