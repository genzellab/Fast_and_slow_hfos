function [signal2]=notch_filter(signal,fs)

%  clear variables 
% fs=1000;
% signal=CORTEX;

%Low pass filter
% Wn=[ 300/(fs/2) ]; % Cutoff=fs_new/2 Hz. 
% [b,a] = butter(3,Wn); %Filter coefficients for LPF. 
% signal=filtfilt(b,a,signal);

%% Compute spectrogram
% [px,f]=prepare_spec(signal,fs);
 
%% Plot spectrogram
%  s=semilogy(f,(px),'Color',[0 0 0],'LineWidth',2);
% 
%  xlabel('Frequency','FontSize',16)
%  ylabel('Log Power','FontSize',16)
%  xlim([0 300])

%% matlab notch filter
n=10; %Filter order
d = designfilt('bandstopiir','FilterOrder',n, ...
               'HalfPowerFrequency1',48,'HalfPowerFrequency2',52, ...
               'DesignMethod','butter','SampleRate',fs);

signal2 = filtfilt(d,signal);

d = designfilt('bandstopiir','FilterOrder',n, ...
               'HalfPowerFrequency1',98,'HalfPowerFrequency2',102, ...
               'DesignMethod','butter','SampleRate',fs);

signal2 = filtfilt(d,signal2);

d = designfilt('bandstopiir','FilterOrder',n, ...
               'HalfPowerFrequency1',148,'HalfPowerFrequency2',152, ...
               'DesignMethod','butter','SampleRate',fs);

signal2 = filtfilt(d,signal2);
% 
d = designfilt('bandstopiir','FilterOrder',n, ...
               'HalfPowerFrequency1',198,'HalfPowerFrequency2',202, ...
               'DesignMethod','butter','SampleRate',fs);

signal2 = filtfilt(d,signal2);
% 
d = designfilt('bandstopiir','FilterOrder',n, ...
               'HalfPowerFrequency1',248,'HalfPowerFrequency2',252, ...
               'DesignMethod','butter','SampleRate',fs);

signal2 = filtfilt(d,signal2);
% 
d = designfilt('bandstopiir','FilterOrder',n, ...
               'HalfPowerFrequency1',298,'HalfPowerFrequency2',302, ...
               'DesignMethod','butter','SampleRate',fs);

           
signal2 = filtfilt(d,signal2);







end