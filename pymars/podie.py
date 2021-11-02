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
#           MARS melibea and CTLearn files.
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
from astropy.time import Time
import uproot


def main():

    parser = argparse.ArgumentParser(
        description=("podie - Python odie-like program to generate the theta2 plot automatically from MARS melibea and CTLearn files."))
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
    parser.add_argument('--times', '-t',
                        help='quate excluded time file; ascii format ".times"')
    args = parser.parse_args()


    # Handling the time slices from Quate to excluded unwanted data
    # TODO_TM: Implement the quate cuts  
    if args.times:
        quate_file = args.times if os.path.isfile(args.times) else args.times.replace(".times","_1.times")
        quate_times = open(quate_file).readlines()
        excluded_time_start, excluded_time_stop = [], []
        for time in quate_times:
            day = time.split()[0].split(".")
            hms = time.split()[1]
            excluded_time_start.append(f"{day[2]}-{day[1]}-{day[0]}T{hms}")
            day = time.split()[2].split(".")
            hms = time.split()[3].replace("\n", "")
            excluded_time_stop.append(f"{day[2]}-{day[1]}-{day[0]}T{hms}")
        excluded_time_start = Time(excluded_time_start, format='isot', scale='utc').mjd
        excluded_time_stop = Time(excluded_time_stop, format='isot', scale='utc').mjd

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
    n_off_per_offregion = np.zeros(len(offRegions))
    
    for i, file in enumerate(files):
        if filename_type == "h5":
            data = pd.HDFStore(file, 'r')
            reco_energy = np.array(data['data']['reco_energy']) * u.TeV
            mc_alt = np.array(data['data']['mc_alt']) * u.rad
            reco_alt = np.array(data['data']['reco_alt']) * u.rad
            mc_az = np.array(data['data']['mc_az']) * u.rad
            reco_az = np.array(data['data']['reco_az']) * u.rad
            pointing_alt = np.array(data['data']['pointing_alt']) * u.rad
            pointing_az = np.array(data['data']['pointing_az']) * u.rad
            print(len(reco_energy))

            mjd = np.array(data['data']['mjd'])
            millisec = np.array(data['data']['millisec'])
            nanosec = np.array(data['data']['nanosec'])

            hadronness = np.array(data['data']['hadronness'])
            size_m1 =  np.array(data['data']['M1_hillas_intensity'])
            size_m2 =  np.array(data['data']['M2_hillas_intensity'])

        elif filename_type == "root":
            melibea_file = uproot.open(file)
            evts = melibea_file["Events"]
            reco_alt = np.asarray(evts["MStereoParDisp.fDirectionY"].array()) +  90.0 - np.asarray(evts["MPointingPos_1.fZd"].array())
            marsDefaultMask = ~np.isnan(reco_alt)
            reco_alt = reco_alt[marsDefaultMask]
            reco_az = np.asarray(evts["MStereoParDisp.fDirectionX"].array())[marsDefaultMask] +  np.asarray(evts["MPointingPos_1.fAz"].array())[marsDefaultMask]
            reco_alt *= u.deg
            reco_az *= u.deg
            mc_alt =  (np.asarray(evts["MSrcPosCam_1.fY"].array())[marsDefaultMask] * 0.00337 +  90.0 - np.asarray(evts["MPointingPos_1.fZd"].array())[marsDefaultMask]) * u.deg
            mc_az =  (np.asarray(evts["MSrcPosCam_1.fX"].array())[marsDefaultMask] * 0.00337 +  np.asarray(evts["MPointingPos_1.fAz"].array())[marsDefaultMask]) * u.deg
            pointing_alt = (90.0 - np.asarray(evts["MPointingPos_1.fZd"].array())[marsDefaultMask]) * u.deg
            pointing_az = np.asarray(evts["MPointingPos_1.fAz"].array())[marsDefaultMask] * u.deg
            reco_energy = np.asarray(evts["MEnergyEst.fEnergy"].array())[marsDefaultMask] * u.GeV

            mjd = np.asarray(evts["MTime_1.fMjd"].array())[marsDefaultMask]
            millisec = np.asarray(evts["MTime_1.fTime.fMilliSec"].array())[marsDefaultMask]/1000.0/3600.0/24.0
            nanosec = np.asarray(evts["MTime_1.fNanoSec"].array())[marsDefaultMask]/1e9/3600.0/24.0

            hadronness = np.asarray(evts["MHadronness.fHadronness"].array())[marsDefaultMask]
            size_m1 = np.asarray(evts["MHillas_1.fSize"].array())[marsDefaultMask]
            size_m2 = np.asarray(evts["MHillas_2.fSize"].array())[marsDefaultMask]

        hadronnessMask = (hadronness < eRange["hadronnessCut"])
        reco_energy = reco_energy[hadronnessMask]
        mc_alt = mc_alt[hadronnessMask]
        reco_alt = reco_alt[hadronnessMask]
        mc_az = mc_az[hadronnessMask]
        reco_az = reco_az[hadronnessMask]
        pointing_alt = pointing_alt[hadronnessMask]
        pointing_az = pointing_az[hadronnessMask]

        mjd = mjd[hadronnessMask]
        millisec = millisec[hadronnessMask]
        nanosec = nanosec[hadronnessMask]

        size_m1 = size_m1[hadronnessMask]
        size_m2 = size_m2[hadronnessMask]

        if eRange["sizeCut"][0] > 0.0 and eRange["sizeCut"][1] > 0.0:
            #TODO: Make generic for any given number of telescopes
            sizeMask = (size_m1 > eRange["sizeCut"][0]) & (size_m2 > eRange["sizeCut"][1])
            reco_energy = reco_energy[sizeMask]
            mc_alt = mc_alt[sizeMask]
            reco_alt = reco_alt[sizeMask]
            mc_az = mc_az[sizeMask]
            reco_az = reco_az[sizeMask]
            pointing_alt = pointing_alt[sizeMask]
            pointing_az = pointing_az[sizeMask]

            mjd = mjd[sizeMask]
            millisec = millisec[sizeMask]
            nanosec = nanosec[sizeMask]

        if eRange["energyCut"] > 0.0:
            energyMask = (reco_energy > eRange["energyCut"])
            reco_energy = reco_energy[energyMask]
            mc_alt = mc_alt[energyMask]
            reco_alt = reco_alt[energyMask]
            mc_az = mc_az[energyMask]
            reco_az = reco_az[energyMask]
            pointing_alt = pointing_alt[energyMask]
            pointing_az = pointing_az[energyMask]

            mjd = mjd[energyMask]
            millisec = millisec[energyMask]
            nanosec = nanosec[energyMask]

        print(len(reco_energy))

        # ON region
        theta2_on = (mc_alt-reco_alt).to(u.deg)**2 + (mc_az-reco_az).to(u.deg)**2
        if i == 0:
            total_theta2_on = theta2_on.value
        else:
            total_theta2_on = np.concatenate((total_theta2_on, theta2_on.value))

        # OFF regions
        delta = np.vstack((mc_alt.to(u.rad)-pointing_alt.to(u.rad), mc_az.to(u.rad)-pointing_az.to(u.rad)))
        for r, rotation in enumerate(offRegions):
            rotation = float(rotation) * u.deg
            rotation_matrix = np.matrix([[np.cos(rotation.to(u.rad)), -np.sin(rotation.to(u.rad))],
                                         [np.sin(rotation.to(u.rad)), np.cos(rotation.to(u.rad))]], dtype=float)
            off = np.squeeze(np.asarray(np.dot(rotation_matrix, delta)))
            off_alt = off[:][0] * u.rad + pointing_alt.to(u.rad)
            off_az = off[:][1] * u.rad + pointing_az.to(u.rad)
            theta2_off = (off_alt-reco_alt).to(u.deg)**2 + (off_az-reco_az).to(u.deg)**2
            n_off_per_offregion[r] += len(theta2_off[theta2_off.value < eRange["signalCut"]])
            if i == 0 and r == 0:
                total_theta2_off = theta2_off.value
            else:
                total_theta2_off = np.concatenate((total_theta2_off, theta2_off.value))
                
    # Create the theta2 plot
    fig, ax = plt.subplots(figsize=(12, 9))

    # Plot the histogram for the ON and OFF region
    ax.hist(total_theta2_on, bins=hist_bin, range=(0.0, args.rangeTh2), color= "black", alpha=1.0, histtype='step', label='On region')
    ax.hist(total_theta2_off, bins=hist_bin, range=(0.0, args.rangeTh2), weights=np.ones(len(total_theta2_off))/float(len(offRegions)), alpha=0.8, label='Off regions', color="grey")

    # Do this from MJD later!
    time = 2.93

    n_on = float(len(total_theta2_on[total_theta2_on < eRange["signalCut"]]))
    total_n_off = np.sum(n_off_per_offregion)
    n_off = total_n_off/float(len(offRegions))
    n_off_error = np.sqrt(n_off/float(len(offRegions)))
    n_excess = n_on - n_off
    n_excess_error = np.sqrt(n_on + n_off_error*n_off_error)
    n_signal = n_on + total_n_off
    alpha = 1.0/float(len(offRegions))
    log_on = np.log( (n_on*(alpha + 1)) / (n_signal*alpha))
    log_off = np.log( (total_n_off*(alpha + 1)) / n_signal)
    sig = np.sqrt((n_on*log_on + total_n_off*log_off) * 2)
    pSigma50 = n_excess / np.sqrt(n_off) * np.sqrt(50.0/time)
    sens = 5./ pSigma50 *100.0
    sens_error = sens * (np.sqrt((1.0/2.0/n_off+1.0/n_excess)*n_off_error*(1.0/2.0/n_off+1.0/n_excess)*n_off_error + n_on/n_excess/n_excess))
    gamma_rate = n_excess/(time*60)
    gamma_rate_error = n_excess_error/(time*60)
    bkg_rate = n_off/(time*60)
    bkg_rate_error = n_off_error/(time*60)

    ax.text(0.8, 0.88, f'Source = Crab Nebula', fontsize=12, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
    ax.text(0.8, 0.84, f'Time = 2.93 h', fontsize=12, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
    ax.text(0.8, 0.8, f'N_on = {int(n_on)}; N_off = {n_off:.1f}\u00B1{n_off_error:.1f} ', fontsize=12, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
    ax.text(0.8, 0.76, f'N_ex = {n_excess:.1f}\u00B1{n_excess_error:.1f} ', fontsize=12, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
    ax.text(0.8, 0.72, f'Significance (Li&Ma) = {sig:.1f}$ \sigma$', fontsize=12, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
    ax.text(0.8, 0.68, f'Sensitivity = {sens:.2f}\u00B1{sens_error:.2f} % Crab', fontsize=12, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
    ax.text(0.8, 0.64, f'Gamma Rate = {gamma_rate:.2f}\u00B1{gamma_rate_error:.2f} / min', fontsize=12, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
    ax.text(0.8, 0.6, f'Bkg Rate = {bkg_rate:.3f}\u00B1{bkg_rate_error:.3f} / min', fontsize=12, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

    ax.set_xbound(0.0, args.rangeTh2)
    ax.axvline(x=eRange["signalCut"], linestyle="dashed", color="red")
    ax.legend(loc='upper center',fontsize=20)

    ax.text(0.5, 0.5, 'Preliminary', fontsize=58, alpha=0.5, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

    ax.set_title(f"{args.name}",fontsize=36)
    ax.set_xlabel(r'$\theta^{2} [deg^{2}]$', fontsize=30)
    ax.set_ylabel(r'$N_{events}$', fontsize=30)
    ax.tick_params(labelsize=25)

    plt.tight_layout()
    # Save plot
    plt.savefig(f"{args.output_dir}/podie_{args.name}.pdf", dpi=600)
    # Show plot
    #plt.show()

if __name__ == "__main__":
    main()


