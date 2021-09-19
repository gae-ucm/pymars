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
#   For help: $ python plotROC.py -h
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

def main():

    parser = argparse.ArgumentParser(
        description=("plotROC - Python ctaplot-based program to plot ROC curves from MARS melibea (under development) and CTLearn files."))
    parser.add_argument('--input_dir', '-i',
                        help='input directory',
                        default="./")
    parser.add_argument('--pattern', '-p',
                        help='pattern to mask unwanted files',
                        default="*")
    parser.add_argument('--output_dir', '-o',
                        help='path where to save generated files. By default, the current directory is used.',
                        default="./")

    args = parser.parse_args()

    # Input handling
    abs_file_dir = os.path.abspath(args.input_dir)
    files = glob.glob(os.path.join(abs_file_dir, args.pattern))
    
    plt.figure(figsize=(12,12))
    ax = plt.axes()
    
    for file in files:
        if file.split('.')[-1] == "h5":

            data = pd.HDFStore(file, 'r')
            for bkg in ['proton', 'off']:
                bkg_mc_particle = np.array(data[bkg]['mc_particle']).astype(int)
                bkg_reco_particle = np.array(data[bkg]['reco_particle']).astype(int)
                bkg_reco_gammaness = np.array(data[bkg]['reco_gammaness']).astype(float)

                gamma_mc_particle = np.array(data['gamma']['mc_particle']).astype(int)
                gamma_reco_particle = np.array(data['gamma']['reco_particle']).astype(int)
                gamma_reco_gammaness = np.array(data['gamma']['reco_gammaness']).astype(float)
                
                # Concatenate the background with the gammas (ringwobble)
                mc_particle = np.concatenate((bkg_mc_particle, gamma_mc_particle))
                reco_particle = np.concatenate((bkg_reco_particle, gamma_reco_particle))
                reco_gammaness = np.concatenate((bkg_reco_gammaness, gamma_reco_gammaness))
                
                fpr, tpr, thresholds = metrics.roc_curve(mc_particle, reco_gammaness)
                            
                ax = ctaplot.plot_roc_curve(mc_particle, reco_gammaness, label="{} (vs {}), AUC={:0.3f}".format(file.split('/')[-1], bkg, metrics.auc(fpr, tpr)), linewidth=1.5)

    ax.set_xlabel('False positive rate', fontsize=25)
    plt.xlim(0.0,1.0)
    ax.set_xticks([0.0,0.2,0.4,0.6,0.8,1.0])

    ax.set_ylabel('True positive rate', fontsize=25)
    plt.ylim(0.0,1.0)
    ax.set_yticks([0.0,0.2,0.4,0.6,0.8,1.0])

    ax.plot([0, 1], [0, 1], '--', color='black')
    #major grid lines
    ax.tick_params(labelsize=25)
    plt.grid(b=True, which='major', color='gray', alpha=0.3, linestyle='dashdot', lw=1.5)

    ax.legend(loc='lower right',fontsize=20)
    ax.set_title("ROC curves",fontsize=30)
    
    plt.savefig(f"{args.output_dir}/roc_curves.png")

if __name__ == "__main__":
    main()
