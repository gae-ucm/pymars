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
#!   Author(s): Tjark Miener, 11/2021 <mailto:tmiener@ucm.es>
#!   Author(s): Alvaro Mas, 11/2021 <mailto:alvmas@ucm.es>
#!
#!   Copyright: pyMARS Software Development, 2021
#!
#!
#\* ================================================================================ */

#######################################################################################
#
#   pulsar - Python pulsar program to generate the pulsar phase plot automatically
#            from MARS melibea and CTLearn files.
#
#   For help: $ pymars-pulsar -h
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
import pint


def main():

    parser = argparse.ArgumentParser(
        description=("pulsar - Python pulsar program to generate the pulsar phase plot automatically from MARS melibea and CTLearn files."))
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

            mjd = np.array(data['data']['mjd'])
            millisec = np.array(data['data']['millisec'])
            nanosec = np.array(data['data']['nanosec'])

            hadronness = np.array(data['data']['hadronness'])
            size_m1 =  np.array(data['data']['M1_hillas_intensity'])
            size_m2 =  np.array(data['data']['M2_hillas_intensity'])
            leakage1_m1 = np.array(data['data']['M1_leakage_intensity_1'])
            leakage1_m2 = np.array(data['data']['M2_leakage_intensity_1'])

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
            leakage1_m1 = np.asarray(evts["MNewImagePar_1.fLeakage1"].array())[marsDefaultMask]
            leakage1_m2 = np.asarray(evts["MNewImagePar_2.fLeakage1"].array())[marsDefaultMask]

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
        leakage1_m1 = leakage1_m1[hadronnessMask]
        leakage1_m2 = leakage1_m2[hadronnessMask]

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

            leakage1_m1 = leakage1_m1[sizeMask]
            leakage1_m2 = leakage1_m2[sizeMask]

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

        # ON region
        theta2_on = (mc_alt-reco_alt).to(u.deg)**2 + (mc_az-reco_az).to(u.deg)**2
        if i == 0:
            total_theta2_on = theta2_on.value
        else:
            total_theta2_on = np.concatenate((total_theta2_on, theta2_on.value))

    # Do this from MJD later!
    #total_time = 2.93

    # Get the events in the ON region 
    on_region = total_theta2_on < eRange["signalCut"]

    mjd = mjd[on_region]
    millisec = millisec[on_region]
    nanosec = nanosec[on_region]

    # From here you have the time_stamps of the events in the ON region:



    # Create the pulsar phase plot
    #fig, ax = plt.subplots(figsize=(12, 9))
    #plt.tight_layout()
    # Save plot
    #plt.savefig(f"{args.output_dir}/pulsar_{args.name}.pdf", dpi=600)
    # Show plot
    #plt.show()

if __name__ == "__main__":
    main()


