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
#   plotRES - Python ctaplot-based program to plot resolutions curves
#             from MARS melibea and CTLearn files.
#
#   For help: $ pymars-plotRES -h
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
        description=("plotRES - Python ctaplot-based program to plot resolutions curves from MARS melibea and CTLearn files."))
    parser.add_argument('--input_dir', '-i',
                        help='input directory',
                        default="./")
    parser.add_argument('--pattern', '-p',
                        help='pattern to mask unwanted files',
                        default="*")
    parser.add_argument('--output_dir', '-o',
                        help='path where to save generated files. By default, the current directory is used.',
                        default="./")
    parser.add_argument('--signalCut', '-s',
                        help='maximum theta2 values for the signal region',
                        default=0.03,
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
    parser.add_argument('--angular', '-a',
                        dest='angular',
                        help='plot angular resolution',
                        default=False,
                        action='store_true')
    parser.add_argument('--energy', '-e',
                        dest='energy',
                        help='plot energy resolution',
                        default=False,
                        action='store_true')
    parser.add_argument('--aleksic',
                        dest='aleksic',
                        help='plot MAGIC performance of Aleksic et al. (2015)',
                        default=False,
                        action='store_true')

    args = parser.parse_args()

    # Input handling
    abs_file_dir = os.path.abspath(args.input_dir)
    files = glob.glob(os.path.join(abs_file_dir, args.pattern))
    
    if args.energy:
        fig, ax = plt.subplots(1,figsize=(17,9))
        for file in files:
            print(file)
            if file.split('.')[-1] == "h5":
                data = pd.HDFStore(file, 'r')
                mc_energy = np.array(data['gamma']['mc_energy']) * u.TeV
                reco_energy = np.array(data['gamma']['reco_energy']) * u.TeV
                mc_alt = np.array(data['gamma']['mc_altitude']) * u.rad
                reco_alt = np.array(data['gamma']['reco_altitude']) * u.rad
                mc_az = np.array(data['gamma']['mc_azimuth']) * u.rad
                reco_az = np.array(data['gamma']['reco_azimuth']) * u.rad

                print(len(reco_energy))

                if args.signalCut or args.hadronnessCut or args.sizeCut or args.leakage1Cut:
                    parameters = np.dstack(data['gamma']["parameters"].to_numpy())
                    hadronness = 1 - data['gamma']["reco_gammaness"].to_numpy()
                    #TODO: Make generic for any given number of telescopes
                    size_m1 = parameters[0][0]
                    size_m2 = parameters[1][0]
                    leakage1_m1 = parameters[0][-2]
                    leakage1_m2 = parameters[1][-2]

                    if args.sizeCut:
                        #TODO: Make generic for any given number of telescopes
                        sizeMask = (size_m1 > args.sizeCut[0]) & (size_m2 > args.sizeCut[1])
                        mc_energy = mc_energy[sizeMask]
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
                        mc_energy = mc_energy[leakage1Mask]
                        reco_energy = reco_energy[leakage1Mask]
                        mc_alt = mc_alt[leakage1Mask]
                        reco_alt = reco_alt[leakage1Mask]
                        mc_az = mc_az[leakage1Mask]
                        reco_az = reco_az[leakage1Mask]
                        hadronness = hadronness[leakage1Mask]

                    if args.signalCut:
                        theta2 = (mc_alt-reco_alt).to(u.deg)**2 + (mc_az-reco_az).to(u.deg)**2
                        theta2Mask = (theta2.value < args.signalCut)
                        mc_energy = mc_energy[theta2Mask]
                        reco_energy = reco_energy[theta2Mask]
                        mc_alt = mc_alt[theta2Mask]
                        reco_alt = reco_alt[theta2Mask]
                        mc_az = mc_az[theta2Mask]
                        reco_az = reco_az[theta2Mask]
                        hadronness = hadronness[theta2Mask]

                    if args.hadronnessCut:
                        hadronnessMask = (hadronness < args.hadronnessCut)
                        mc_energy = mc_energy[hadronnessMask]
                        reco_energy = reco_energy[hadronnessMask]
                        mc_alt = mc_alt[hadronnessMask]
                        reco_alt = reco_alt[hadronnessMask]
                        mc_az = mc_az[hadronnessMask]
                        reco_az = reco_az[hadronnessMask]

                    print(len(reco_energy))
            elif file.split('.')[-1] == "root":
                melibea_file = uproot.open(file)
                evts = melibea_file["Events"]
                reco_alt = np.asarray(evts["MStereoParDisp.fDirectionY"].array()) +  90.0 - np.asarray(evts["MPointingPos_1.fZd"].array())
                reco_az = np.asarray(evts["MStereoParDisp.fDirectionX"].array()) +  np.asarray(evts["MPointingPos_1.fAz"].array())
                marsDefaultMask = ~np.isnan(reco_alt)
                reco_alt *= u.deg
                reco_az *= u.deg
                mc_alt =  (np.asarray(evts["MSrcPosCam_1.fY"].array()) * 0.00337 +  90.0 - np.asarray(evts["MPointingPos_1.fZd"].array())) * u.deg
                mc_az =  (np.asarray(evts["MSrcPosCam_1.fX"].array()) * 0.00337 +  np.asarray(evts["MPointingPos_1.fAz"].array())) * u.deg
                mc_energy =  np.asarray(evts["MMcEvt_1.fEnergy"].array()) * u.GeV
                reco_energy = np.asarray(evts["MEnergyEst.fEnergy"].array()) * u.GeV

                print(len(reco_energy))

                if args.signalCut or args.hadronnessCut or args.sizeCut or args.leakage1Cut:
                    mc_energy = mc_energy[marsDefaultMask]
                    reco_energy = reco_energy[marsDefaultMask]
                    mc_alt = mc_alt[marsDefaultMask]
                    reco_alt = reco_alt[marsDefaultMask]
                    mc_az = mc_az[marsDefaultMask]
                    reco_az = reco_az[marsDefaultMask]
                    hadronness = np.asarray(evts["MHadronness.fHadronness"].array())[marsDefaultMask]
                    size_m1 = np.asarray(evts["MHillas_1.fSize"].array())[marsDefaultMask]
                    size_m2 = np.asarray(evts["MHillas_2.fSize"].array())[marsDefaultMask]
                    leakage1_m1 = np.asarray(evts["MNewImagePar_1.fLeakage1"].array())[marsDefaultMask]
                    leakage1_m2 = np.asarray(evts["MNewImagePar_2.fLeakage1"].array())[marsDefaultMask]

                    if args.sizeCut:
                        sizeMask = (size_m1 > args.sizeCut[0]) & (size_m2 > args.sizeCut[1])
                        mc_energy = mc_energy[sizeMask]
                        reco_energy = reco_energy[sizeMask]
                        mc_alt = mc_alt[sizeMask]
                        reco_alt = reco_alt[sizeMask]
                        mc_az = mc_az[sizeMask]
                        reco_az = reco_az[sizeMask]
                        hadronness = hadronness[sizeMask]
                        leakage1_m1 = leakage1_m1[sizeMask]
                        leakage1_m2 = leakage1_m2[sizeMask]

                    if args.leakage1Cut:
                        leakage1Mask = (leakage1_m1 < args.leakage1Cut[0]) & (leakage1_m2 < args.leakage1Cut[1])
                        mc_energy = mc_energy[leakage1Mask]
                        reco_energy = reco_energy[leakage1Mask]
                        mc_alt = mc_alt[leakage1Mask]
                        reco_alt = reco_alt[leakage1Mask]
                        mc_az = mc_az[leakage1Mask]
                        reco_az = reco_az[leakage1Mask]
                        hadronness = hadronness[leakage1Mask]

                    if args.signalCut:
                        theta2 = (mc_alt-reco_alt).to(u.deg)**2 + (mc_az-reco_az).to(u.deg)**2
                        theta2Mask = (theta2.value < args.signalCut)
                        mc_energy = mc_energy[theta2Mask]
                        reco_energy = reco_energy[theta2Mask]
                        mc_alt = mc_alt[theta2Mask]
                        reco_alt = reco_alt[theta2Mask]
                        mc_az = mc_az[theta2Mask]
                        reco_az = reco_az[theta2Mask]
                        hadronness = hadronness[theta2Mask]

                    if args.hadronnessCut:
                        hadronnessMask = (hadronness < args.hadronnessCut)
                        mc_energy = mc_energy[hadronnessMask]
                        reco_energy = reco_energy[hadronnessMask]
                        mc_alt = mc_alt[hadronnessMask]
                        reco_alt = reco_alt[hadronnessMask]
                        mc_az = mc_az[hadronnessMask]
                        reco_az = reco_az[hadronnessMask]

                    print(len(reco_energy))

            ax = ctaplot.plot_energy_resolution(mc_energy.to(u.TeV).value, reco_energy.to(u.TeV).value, label=f"{file.split('/')[-1]}")            

        if args.aleksic:
            # Aleksic et al. (2015) MAGIC performance (energy)
            aleksic_energy = np.array([60.14412430231677, 88.92924106247686, 146.42158223581913, 378.0025174064791, 598.5081431652516, 1500.4475087782453, 2329.7193634991513, 3761.5908025786603, 5955.893474058259, 9430.230170174897, 14786.02308938286])*1e-3
            aleksic_resolution = np.array([0.21699999999999994, 0.20099999999999998, 0.18, 0.15499999999999997, 0.14799999999999996, 0.16099999999999998, 0.18, 0.19599999999999998, 0.21899999999999997, 0.22699999999999995, 0.20799999999999996])
            ax.plot(aleksic_energy, aleksic_resolution, label="Aleksic et al. (2015)")

        # Scale, labels and title
        ax.set_ylabel(r"$(\Delta E/E)_{68}$",fontsize=25)

        ax.set_ybound(0.08,0.33)
        ax.set_yticks([0.1,0.15,0.2,0.25,0.3])

        ax.set_xlabel("Energy [TeV]",fontsize=25)
        ax.set_xscale('log')
        ax.set_xbound(4e-2,25)
        ax.set_title('Energy resolution',fontsize=30)

        ax.tick_params(labelsize=25)

        #major and minor grid lines
        plt.grid(b=True, which='major', color='gray', alpha=0.3, linestyle='dashdot', lw=1.5)
        plt.minorticks_on()
        plt.grid(b=True, which='minor', color='beige', alpha=0.5, ls='-', lw=1)
        ax.legend(loc='upper right',fontsize=25)

        plt.savefig(f"{args.output_dir}/energyresolution.png")
    
    if args.angular:
        fig, ax = plt.subplots(1,figsize=(17,9))
    
        for file in files:
            print(file)
            if file.split('.')[-1] == "h5":
                data = pd.HDFStore(file, 'r')
                reco_energy = np.array(data['gamma']['reco_energy']) * u.TeV
                mc_alt = np.array(data['gamma']['mc_altitude']) * u.rad
                reco_alt = np.array(data['gamma']['reco_altitude']) * u.rad
                mc_az = np.array(data['gamma']['mc_azimuth']) * u.rad
                reco_az = np.array(data['gamma']['reco_azimuth']) * u.rad
                #if args.signalCut or args.hadronnessCut or args.sizeCut or args.leakage1Cut:
                print(len(reco_energy))

                if args.hadronnessCut or args.sizeCut or args.leakage1Cut:
                    parameters = np.dstack(data['gamma']["parameters"].to_numpy())
                    hadronness = 1 - data['gamma']["reco_gammaness"].to_numpy()
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

                    #if args.signalCut:
                    #    theta2 = (mc_alt-reco_alt).to(u.deg)**2 + (mc_az-reco_az).to(u.deg)**2
                    #    theta2Mask = (theta2.value < args.signalCut)
                    #    mc_energy = mc_energy[theta2Mask]
                    #    reco_energy = reco_energy[theta2Mask]
                    #    mc_alt = mc_alt[theta2Mask]
                    #    reco_alt = reco_alt[theta2Mask]
                    #    mc_az = mc_az[theta2Mask]
                    #    reco_az = reco_az[theta2Mask]
                    #    hadronness = hadronness[theta2Mask]

                    if args.hadronnessCut:
                        hadronnessMask = (hadronness < args.hadronnessCut)
                        reco_energy = reco_energy[hadronnessMask]
                        mc_alt = mc_alt[hadronnessMask]
                        reco_alt = reco_alt[hadronnessMask]
                        mc_az = mc_az[hadronnessMask]
                        reco_az = reco_az[hadronnessMask]

                    print(len(reco_energy))

            elif file.split('.')[-1] == "root":
                melibea_file = uproot.open(file)
                evts = melibea_file["Events"]
                reco_alt = np.asarray(evts["MStereoParDisp.fDirectionY"].array()) +  90.0 - np.asarray(evts["MPointingPos_1.fZd"].array())
                reco_az = np.asarray(evts["MStereoParDisp.fDirectionX"].array()) +  np.asarray(evts["MPointingPos_1.fAz"].array())
                marsDefaultMask = ~np.isnan(reco_alt)
                reco_alt *= u.deg
                reco_az *= u.deg
                mc_alt =  (np.asarray(evts["MSrcPosCam_1.fY"].array()) * 0.00337 +  90.0 - np.asarray(evts["MPointingPos_1.fZd"].array())) * u.deg
                mc_az =  (np.asarray(evts["MSrcPosCam_1.fX"].array()) * 0.00337 +  np.asarray(evts["MPointingPos_1.fAz"].array())) * u.deg
                reco_energy = np.asarray(evts["MEnergyEst.fEnergy"].array()) * u.GeV

                print(len(reco_energy))
                #if args.signalCut or args.hadronnessCut or args.sizeCut or args.leakage1Cut:

                if args.hadronnessCut or args.sizeCut or args.leakage1Cut:
                    reco_energy = reco_energy[marsDefaultMask]
                    mc_alt = mc_alt[marsDefaultMask]
                    reco_alt = reco_alt[marsDefaultMask]
                    mc_az = mc_az[marsDefaultMask]
                    reco_az = reco_az[marsDefaultMask]
                    hadronness = np.asarray(evts["MHadronness.fHadronness"].array())[marsDefaultMask]
                    size_m1 = np.asarray(evts["MHillas_1.fSize"].array())[marsDefaultMask]
                    size_m2 = np.asarray(evts["MHillas_2.fSize"].array())[marsDefaultMask]
                    leakage1_m1 = np.asarray(evts["MNewImagePar_1.fLeakage1"].array())[marsDefaultMask]
                    leakage1_m2 = np.asarray(evts["MNewImagePar_2.fLeakage1"].array())[marsDefaultMask]

                    if args.sizeCut:
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
                        leakage1Mask = (leakage1_m1 < args.leakage1Cut[0]) & (leakage1_m2 < args.leakage1Cut[1])
                        reco_energy = reco_energy[leakage1Mask]
                        mc_alt = mc_alt[leakage1Mask]
                        reco_alt = reco_alt[leakage1Mask]
                        mc_az = mc_az[leakage1Mask]
                        reco_az = reco_az[leakage1Mask]
                        hadronness = hadronness[leakage1Mask]

                    #if args.signalCut:
                    #    theta2 = (mc_alt-reco_alt).to(u.deg)**2 + (mc_az-reco_az).to(u.deg)**2
                    #    theta2Mask = (theta2.value < args.signalCut)
                    #    reco_energy = reco_energy[theta2Mask]
                    #    mc_alt = mc_alt[theta2Mask]
                    #    reco_alt = reco_alt[theta2Mask]
                    #    mc_az = mc_az[theta2Mask]
                    #    reco_az = reco_az[theta2Mask]
                    #    hadronness = hadronness[theta2Mask]

                    if args.hadronnessCut:
                        hadronnessMask = (hadronness < args.hadronnessCut)
                        reco_energy = reco_energy[hadronnessMask]
                        mc_alt = mc_alt[hadronnessMask]
                        reco_alt = reco_alt[hadronnessMask]
                        mc_az = mc_az[hadronnessMask]
                        reco_az = reco_az[hadronnessMask]
                    print(len(reco_energy))

            ax = ctaplot.plot_angular_resolution_per_energy(reco_alt.to(u.rad).value, reco_az.to(u.rad).value, mc_alt.to(u.rad).value, mc_az.to(u.rad).value, reco_energy.to(u.TeV).value, bias_correction=False, label=f"{file.split('/')[-1]}")
            
        if args.aleksic:
            # Aleksic et al. (2015) MAGIC performance (angular)
            aleksic_energy = np.array([80.27961169117373, 98.4965576498376, 149.39707865342058, 235.34800775765504, 379.26901907322497, 597.4695679970707, 948.3609962468593, 1493.9707865342073, 2389.4033799062745, 3792.6901907322535, 6020.121593466931, 9555.714328759786])*1e-3
            aleksic_resolution = np.array([0.16343127141077546, 0.15305784802242292, 0.12609730224229213, 0.10244695966988478, 0.08252238399252569, 0.06881286982248522, 0.0613156921519776, 0.05547600046714421, 0.05419115929616944, 0.05249279819370914, 0.05079443709124884, 0.054480691373403944])
            ax.plot(aleksic_energy, aleksic_resolution, label="Aleksic et al. (2015)")


        # Scale, labels and title
        ax.set_ylabel(r'$\theta [deg]$',fontsize=25)

        ax.set_ybound(0.03,0.22)
        ax.set_yticks([0.05,0.1,0.15,0.2])

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

        plt.savefig(f"{args.output_dir}/angularresolution.png")

if __name__ == "__main__":
    main()


