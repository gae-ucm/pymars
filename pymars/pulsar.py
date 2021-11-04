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
from pint import toa
from pint import models
import pint
import csv

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
    parser.add_argument('--energyCut', '-e',
                        help='minimum reco energy values',
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
    files = np.sort(glob.glob(os.path.join(abs_file_dir, args.pattern)))
    filename_type = files[0].split('.')[-1]

    phases = np.array([])
    for i, file in enumerate(files):
        print(file)
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
            print(np.asarray(evts["MTime_1.fTime.fMilliSec"].array())[marsDefaultMask])
            hadronness = np.asarray(evts["MHadronness.fHadronness"].array())[marsDefaultMask]
            size_m1 = np.asarray(evts["MHillas_1.fSize"].array())[marsDefaultMask]
            size_m2 = np.asarray(evts["MHillas_2.fSize"].array())[marsDefaultMask]
            leakage1_m1 = np.asarray(evts["MNewImagePar_1.fLeakage1"].array())[marsDefaultMask]
            leakage1_m2 = np.asarray(evts["MNewImagePar_2.fLeakage1"].array())[marsDefaultMask]

        hadronnessMask = (hadronness < args.hadronnessCut)
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

        if args.sizeCut[0] > 0.0 and args.sizeCut[1] > 0.0:
            #TODO: Make generic for any given number of telescopes
            sizeMask = (size_m1 > args.sizeCut[0]) & (size_m2 > args.sizeCut[1])
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

        if args.energyCut > 0.0:
            energyMask = (reco_energy > args.energyCut)
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
        # Do this from MJD later!
        #total_time = 2.93

        # Get the events in the ON region 
        on_region = theta2_on.value < args.signalCut

        mjd = mjd[on_region]
        millisec = millisec[on_region]
        nanosec = nanosec[on_region]
        timestamp = np.float128(mjd+millisec+nanosec) 
 
        # From here you have the time_stamps of the events in the ON region:        
	
	# Write the .tim file to use with PINT
        timname='times.tim'
        timFile=open(timname,'w+')
        timFile.write('FORMAT 1 \n')
        for i in range(0,len(timestamp)):
                timFile.write('magic '+'0.0 '+str(timestamp[i])+' 0.0 '+ 'magic'+' \n')
        timFile.close()

        #Upload the TOAs object
        t= pint.toa.get_TOAs(timname)
    
        #Create model from ephemeris
        #Read the ephemeris txt file
        colnames=['PSR', 'RAJ1','RAJ2','RAJ3', 'DECJ1','DECJ2','DECJ3', 'START', 'FINISH', 't0geo', 'F0', 'F1', 'F2',
              'RMS','Observatory', 'EPHEM', 'PSR2']
        df_ephem = pd.read_csv(abs_file_dir+'/all.gro', delimiter='\s+',names=colnames,header=None)    
        
        #Search the line of the ephemeris at which the interval time of arrivals given belongs
        for i in range(0,len(df_ephem['START'])):
                if (timestamp[0]>df_ephem['START'][i]) & (timestamp[0]<df_ephem['FINISH'][i]):
                        break		
                elif (timestamp[0]< df_ephem['START'][i]) & (i==0):
                        print('No ephemeris available')
                elif (timestamp[0]> df_ephem['START'][i]) & (timestamp[0]> df_ephem['FINISH'][i])& (i==len(df_ephem['START'])):
                        print('No ephemeris available')

       

        # Select the components(see PINT pulsar documentation)
        components=[]
        for name in ["AbsPhase","AstrometryEquatorial", "Spindown","SolarSystemShapiro"]:
                component_object = models.timing_model.Component.component_types[name]  
                components.append(component_object())

        time_model = models.timing_model.TimingModel(components=components)


        #Add second derivative of the frequency
        f2 = models.parameter.prefixParameter(
                parameter_type="float",
                name="F2",
                value=0.0,
                units=u.Hz / (u.s) ** 2,
                longdouble=True,
        )	
        time_model.components["Spindown"].add_param(f2, setup=True)


        #Add START and FINISH parameters
        time_model.add_param_from_top(models.parameter.MJDParameter(name="START", description="Start MJD for fitting"), "")
        time_model.add_param_from_top(models.parameter.MJDParameter(name="FINISH", description="End MJD for fitting"), "")
    
	#Adopt format of f1 and f2 from PINT        
        f1=float(str(df_ephem['F1'][i].replace('D','E')))
        f2=float(str(df_ephem['F2'][i].replace('D','E')))
    
        #Create a dictionary with the values of the parameters
        param_dic = {
                "PSR":(df_ephem['PSR'][i]) ,
                "RAJ": (str(df_ephem['RAJ1'][i])+':'+ str(df_ephem['RAJ2'][i])+':'+str(df_ephem['RAJ3'][i])),
                "DECJ": (str(df_ephem['DECJ1'][i])+':'+ str(df_ephem['DECJ2'][i])+':'+str(df_ephem['DECJ3'][i])),
                "START": (Time(df_ephem['START'][i], format="mjd", scale="tdb")),
                "FINISH": (Time(df_ephem['FINISH'][i], format="mjd", scale="tdb")),
                "EPHEM":(df_ephem['EPHEM'][i]),
       	        'PEPOCH':(Time(int(df_ephem['t0geo'][i]), format="mjd", scale="tdb")),
                "F0": (df_ephem['F0'][i]*u.Hz),
                "F1": (f1*u.Hz/u.s),
                "F2":(f2*u.Hz/(u.s**2)),
                "TZRMJD":(Time(df_ephem['t0geo'][i], format="mjd", scale="tdb")),
                "TZRFRQ":(0.0*u.Hz),
                "TZRSITE":('coe'),
                }

	
        #Create the model using PINT
        for name_par, value in param_dic.items():
                p = getattr(time_model, name_par)
                p.quantity = value
                
        time_model.validate()
  
        print(time_model)

        #Create the .par file
        parname="model.par"
        f=open(parname,"w+")
        f.write(time_model.as_parfile())
        f.close()

        #Upload TOAs and model
        m=models.get_model(parname)

        #Calculate the phases
        print('Calculating barycentric time and absolute phase')
        barycent_toas=m.get_barycentric_toas(t)
        phase_tuple=m.phase(t,abs_phase=True)
        phases=np.concatenate([phases,phase_tuple.frac])

	#Remove the tim and par files
        os.remove(str(os.getcwd())+'/'+timname)
        os.remove(str(os.getcwd())+'/'+parname)

    
    #Shift phases
    for i in range(0,len(phases)):
        if phases[i]<0:
             phases[i]=phases[i]+1


    #Statistics
    limitP1=[0.026,0.983]
    limitP2=[0.377,0.422]
    limitOFF=[0.52,0.87]

    P1_counts=len(phases[(phases<limitP1[0]) | (phases>limitP1[1])])
    P2_counts=len(phases[(phases>limitP2[0]) & (phases<limitP2[1])])    
    OFF_counts=len(phases[(phases>limitOFF[0]) & (phases<limitOFF[1])])

    deltaP1=1+limitP1[0]-limitP1[1]
    deltaP2=limitP2[1]-limitP2[0]
    deltaOFF=limitOFF[1]-limitOFF[0]

    alpha=deltaP1/deltaOFF
    S1=np.sqrt(2)*(P1_counts*np.log((1+alpha)/alpha*(P1_counts/(P1_counts+OFF_counts)))+OFF_counts*np.log((1+alpha)*(OFF_counts/(P1_counts+OFF_counts))))**(1/2)
    alpha=deltaP2/deltaOFF
    S2=np.sqrt(2)*(P2_counts*np.log((1+alpha)/alpha*(P2_counts/(P2_counts+OFF_counts)))+OFF_counts*np.log((1+alpha)*(OFF_counts/(P2_counts+OFF_counts))))**(1/2)


    #Histogram of phases
    lc = np.histogram(phases, np.linspace(0,1,21))    
    pcentres=(lc[1][1:]+lc[1][:-1])/2
   
    # Create the pulsar phase plot
    fig, ax = plt.subplots(figsize=(12, 9))
    plt.bar(pcentres,lc[0],width=1/len(lc[0]),color='blue',alpha=0.5,edgecolor = 'black',linewidth=2)
    plt.bar(pcentres+np.ones(len(pcentres)),lc[0],width=1/len(lc[0]),color='blue',alpha=0.5,edgecolor ='black',linewidth=2)
        
    #Plot errorbars
    plt.errorbar(pcentres,lc[0],yerr=np.sqrt(lc[0]),color='black',fmt='.')
    plt.errorbar(pcentres+np.ones(len(pcentres)),lc[0],yerr=np.sqrt(lc[0]),color='black',fmt='.')
    

    plt.fill_between(np.linspace(0,limitP1[0],150), 0, 1600500,facecolor='orange',color='orange',alpha=0.2)
    plt.fill_between(np.linspace(limitP1[1],1+limitP1[0],150), 0, 1600500,facecolor='orange',color='orange',alpha=0.2)
    plt.fill_between(np.linspace(1+limitP1[1],2,150), 0, 1600500,facecolor='orange',color='orange',alpha=0.2)

    plt.fill_between(np.linspace(limitP2[0],limitP2[1],150), 0, 1600500,facecolor='green',color='green',alpha=0.2)    
    plt.fill_between(np.linspace(limitP2[0]+1,limitP2[1]+1,150), 0, 1600500,facecolor='green',color='green',alpha=0.2)

    #Add labels
    plt.xlabel('Pulsar phase')
    plt.ylabel('Events')
    text_towrite=f'P1:Sig(Li&Ma):{S1:.2f}$\sigma$'+'\n'+f'P2:Sig(Li&Ma):{S2:.2f}$\sigma$'
    #Set limits in axis
    plt.ylim(max(min(lc[0])-3*np.sqrt(min(lc[0])),0),max(lc[0])+2*np.sqrt(max(lc[0])))
    plt.annotate(text_towrite, xy=(0.7, 0.9), xytext=(0.7,0.9), fontsize=15,xycoords='axes fraction', textcoords='offset points', color='tab:red',bbox=dict(facecolor='white', edgecolor='tab:red'),horizontalalignment='left', verticalalignment='top')
        
    plt.tight_layout()
    # Save plot
    plt.savefig(f"{args.output_dir}/pulsar_{args.name}.pdf", dpi=600)
    # Show plot
    plt.show()

if __name__ == "__main__":
    main()





