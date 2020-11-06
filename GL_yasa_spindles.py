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
aver=mat['V_par']

averout=aver
ajalas=np.zeros((1,aver.size))

for x in range(0, aver.size):
    data=aver[x];
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
    
    
    if sp is not None:
        ajalas[0,x]=1
       #  averout[x,0]=sp.get_mask()
        summary=sp.summary()
        ch=pd.DataFrame(summary).to_numpy()

        averout[x,0]=ch
    else:
        averout[x,0]=[]
        

scipy.io.savemat('YASA_PAR_spindles.mat',{'averout':averout})



mat=scipy.io.loadmat('YASA_PFC.mat')
aver=mat['V_pfc']

averout=aver
ajalas=np.zeros((1,aver.size))

for x in range(0, aver.size):
    data=aver[x];
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
        ajalas[0,x]=1
       #  averout[x,0]=sp.get_mask()
        summary=sp.summary()
        ch=pd.DataFrame(summary).to_numpy()

        averout[x,0]=ch
    else:
        averout[x,0]=[]

scipy.io.savemat('YASA_PFC_spindles.mat',{'averout':averout})


 #summary=sp.summary()
 #ch=pd.DataFrame(summary).to_numpy()
