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
#   plotROC - Python ctaplot-based program to plot ROC curves
#             from MARS melibea (under development) and CTLearn files.
#
#   For help: $ pymars-plotROC -h
#
#######################################################################################

import argparse
import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import metrics
import ctaplot
import uproot

def main():

    parser = argparse.ArgumentParser(
        description=("plotROC - Python ctaplot-based program to plot ROC curves from MARS melibea (under development) and CTLearn files."))
    parser.add_argument('--input_dir', '-i',
                        help='input directory',
                        default="./")
    parser.add_argument('--pattern', '-p',
                        help='pattern to mask unwanted files',
                        default="*")
    parser.add_argument('--input_offdata', '-j',
                        help='input directory for the off data (melibea)',
                        default=None)
    parser.add_argument('--pattern_offdata', '-q',
                        help='pattern to mask unwanted files for the off data input directory',
                        default="*")
    parser.add_argument('--output_dir', '-o',
                        help='path where to save generated files. By default, the current directory is used.',
                        default="./")

    args = parser.parse_args()

    # Input handling
    abs_file_dir = os.path.abspath(args.input_dir)
    files = glob.glob(os.path.join(abs_file_dir, args.pattern))
    if args.input_offdata:
        abs_file_dir_off = os.path.abspath(args.input_offdata)
        off_files = glob.glob(os.path.join(abs_file_dir_off, args.pattern_offdata))


    fig, ax = plt.subplots(1,figsize=(12,12))

    # color cycle of matplotlib
    color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
    color_cycle = ["orange", "royalblue"]
    naming = ["CTLearn - cleaned images", "CTLearn - raw images"]

    color_counter = 0    
    for file in files:
        if file.split('.')[-1] == "h5":
            color_counter += 1
            color_counter %= len(color_cycle)
            data = pd.HDFStore(file, 'r')
            linestyle = {'proton': 'dashed', 'off': 'solid'}
            for bkg in ['proton', 'off']:
                bkg_mc_particle = np.array(data[bkg]['mc_particle']).astype(int)
                bkg_reco_gammaness = np.array(data[bkg]['reco_gammaness']).astype(float)

                gamma_mc_particle = np.array(data['gamma']['mc_particle']).astype(int)
                gamma_reco_gammaness = np.array(data['gamma']['reco_gammaness']).astype(float)
                 
                # Concatenate the background with the gammas (ringwobble)
                mc_particle = np.concatenate((bkg_mc_particle, gamma_mc_particle))
                reco_gammaness = np.concatenate((bkg_reco_gammaness, gamma_reco_gammaness))
                
                fpr, tpr, thresholds = metrics.roc_curve(mc_particle, reco_gammaness)
                particle_type = "MCs proton" if bkg == 'proton' else "off data"
                ax.plot(fpr, tpr, label="{} (vs {}), AUC={:0.3f}".format(naming[color_counter], particle_type, metrics.auc(fpr, tpr)), color=color_cycle[color_counter], linestyle=linestyle[bkg], linewidth=1.5)
                #ax = ctaplot.plot_roc_curve(mc_particle, reco_gammaness, label="{} (vs {}), AUC={:0.3f}".format(naming[color_counter], particle_type, metrics.auc(fpr, tpr)), color=color_cycle[color_counter], linestyle=linestyle[bkg], linewidth=1.5)

        elif file.split('.')[-1] == "root":
            melibea_file = uproot.open(file)
            evts = melibea_file["Events"]
            if b'MMcEvt_1.' in evts.keys():
                if int(evts["MMcEvt_1.fPartId"].array()[0]) == 1:
                    gamma_mc_particle_root = np.asarray(evts["MMcEvt_1.fPartId"].array())
                    gamma_reco_gammaness_root = 1.0 - np.asarray(evts["MHadronness.fHadronness"].array())

    for file in files:
        if file.split('.')[-1] == "root":
            melibea_file = uproot.open(file)
            evts = melibea_file["Events"]
            if b'MMcEvt_1.' not in evts.keys(): continue
            if int(evts["MMcEvt_1.fPartId"].array()[0]) == 1: continue

            bkg_reco_gammaness = 1.0 - np.asarray(evts["MHadronness.fHadronness"].array())
            bkg_mc_particle = np.zeros(len(bkg_reco_gammaness))

            # Concatenate the background with the gammas (ringwobble) from MARS
            mc_particle = np.concatenate((bkg_mc_particle, gamma_mc_particle_root))
            reco_gammaness = np.concatenate((bkg_reco_gammaness, gamma_reco_gammaness_root))

            fpr, tpr, thresholds = metrics.roc_curve(mc_particle, reco_gammaness)
            ax.plot(fpr, tpr, label="MARS - standard RF (vs MCs proton), AUC={:0.3f}".format(metrics.auc(fpr, tpr)), color='black', linestyle="dashed", linewidth=1.5)
            #ax = ctaplot.plot_roc_curve(mc_particle, reco_gammaness, label="MARS - standard RF (vs MCs proton), AUC={:0.3f}".format(metrics.auc(fpr, tpr)), color='black', linestyle="dashed", linewidth=1.5)

    if args.input_offdata:
        bkg_mc_particle = []
        bkg_reco_gammaness = []
        for file in off_files:
            if file.split('.')[-1] == "root":
                melibea_file = uproot.open(file)
                evts = melibea_file["Events"]
                if b'MMcEvt_1.' in evts.keys(): continue
                file_reco_gammaness = 1.0 - np.asarray(evts["MHadronness.fHadronness"].array())
                file_mc_particle = np.zeros(len(file_reco_gammaness))
                # Concatenate the background with the gammas (ringwobble)
                bkg_mc_particle = np.concatenate((bkg_mc_particle, file_mc_particle))
                bkg_reco_gammaness = np.concatenate((bkg_reco_gammaness, file_reco_gammaness))

        # Concatenate the background with the gammas (ringwobble) from MARS
        mc_particle = np.concatenate((bkg_mc_particle, gamma_mc_particle_root))
        reco_gammaness = np.concatenate((bkg_reco_gammaness, gamma_reco_gammaness_root))

        fpr, tpr, thresholds = metrics.roc_curve(mc_particle, reco_gammaness)
        ax.plot(fpr, tpr, label="MARS - standard RF (vs off data), AUC={:0.3f}".format(metrics.auc(fpr, tpr)), color='black', linestyle="solid", linewidth=1.5)
        #ax = ctaplot.plot_roc_curve(mc_particle, reco_gammaness, label="MARS - standard RF (vs off data), AUC={:0.3f}".format(metrics.auc(fpr, tpr)), color='black', linestyle="solid", linewidth=1.5)

    ax.set_xlabel('False positive rate', fontsize=34)
    ax.set_xbound(0.0,1.0)
    ax.set_xticks([0.0,1.0])

    ax.set_ylabel('True positive rate', fontsize=34)
    ax.set_ybound(0.0,1.0)
    ax.set_yticks([0.0,1.0])

    ax.plot([0, 1], [0, 1], '--', color='grey')
    #major grid lines
    ax.tick_params(labelsize=30)
    #plt.grid(b=True, which='major', color='gray', alpha=0.3, linestyle='dashdot', lw=1.5)

    ax.text(0.5, 0.5, 'Preliminary', fontsize=58, alpha=0.5, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

    ax.legend(loc='lower right',fontsize=22)
    ax.set_title("ROC curves",fontsize=40)
    plt.tight_layout()
    
    plt.savefig(f"{args.output_dir}/roc_curves.pdf", dpi=600)

if __name__ == "__main__":
    main()
