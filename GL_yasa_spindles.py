# -*- coding: utf-8 -*-

import yasa
import scipy.io
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

sns.set(font_scale=1.2)

#data = np.loadtxt('data_spindles_15s_200Hz.txt')
mat=scipy.io.loadmat('YASA_PAR.mat')
bouts=mat['V_hpc']

boutsout=bouts
about=np.zeros((1,bouts.size))

#Iterate across all bouts
for x in range(0, bouts.size):
    data=bouts[x];
    data=data[0];
    
    data=np.transpose(data)
    data=data[0];
    
    # Define sampling frequency and time vector
    sf = 1000.
    times = np.arange(data.size) / sf
    
    # # Plot the signal
    # fig, ax = plt.subplots(1, 1, figsize=(14, 4))
    # plt.plot(times, data, lw=1.5, color='k')
    # plt.xlabel('Time (seconds)')
    # plt.ylabel('Amplitude (uV)')
    # plt.xlim([times.min(), times.max()])
    # plt.title('LFP NREM epoch')
    # sns.despine()
    
    if data.size>sf :
        sp = yasa.spindles_detect(data, sf);
    else:
        print('Short epoch')
        sp = None;
    
    
    if sp is not None:  #If there was a spindle detected
        about[0,x]=1
       #  boutsout[x,0]=sp.get_mask()
        summary=sp.summary()
        ch=pd.DataFrame(summary).to_numpy()

        boutsout[x,0]=ch
    else:
        boutsout[x,0]=[]
        

scipy.io.savemat('YASA_PAR_spindles.mat',{'boutsout':boutsout})



mat=scipy.io.loadmat('YASA_PFC.mat')
bouts=mat['V_pfc']

boutsout=bouts
about=np.zeros((1,bouts.size))

for x in range(0, bouts.size):
    data=bouts[x];
    data=data[0];
    
    data=np.transpose(data)
    data=data[0];
    
    # Define sampling frequency and time vector
    sf = 1000.
    times = np.arange(data.size) / sf
    
    # # Plot the signal
    # fig, ax = plt.subplots(1, 1, figsize=(14, 4))
    # plt.plot(times, data, lw=1.5, color='k')
    # plt.xlabel('Time (seconds)')
    # plt.ylabel('Amplitude (uV)')
    # plt.xlim([times.min(), times.max()])
    # plt.title('LFP NREM epoch')
    # sns.despine()
    if data.size>sf :
        sp = yasa.spindles_detect(data, sf);
    else:
        print('Short epoch')
        sp = None;
        
    #sp = yasa.spindles_detect(data, sf)
    if sp is not None:
        about[0,x]=1
       #  boutsout[x,0]=sp.get_mask()
        summary=sp.summary()
        ch=pd.DataFrame(summary).to_numpy()

        boutsout[x,0]=ch
    else:
        boutsout[x,0]=[]

scipy.io.savemat('YASA_PFC_spindles.mat',{'boutsout':boutsout})


 #summary=sp.summary()
 #ch=pd.DataFrame(summary).to_numpy()
