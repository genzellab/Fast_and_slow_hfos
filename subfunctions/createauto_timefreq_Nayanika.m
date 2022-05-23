function [granger, freq]=createauto_timefreq_Nayanika(data1,freqrange)
% Compute Multitaper
equis=0.5;
       
cfg = [];
cfg.method = 'mtmconvol';
%cfg.method = 'wavelet';

cfg.taper = 'hanning';
%cfg.taper = 'dpss';

%cfg.pad='nextpow2';
cfg.pad=3;
%cfg.padtype='edge';
cfg.foi = freqrange;
%cfg.t_ftimwin = 0.4 * ones(size(cfg.foi)); 
cfg.t_ftimwin = 0.5.*ones(size(4./cfg.foi')); 
cfg.tapsmofrq = equis*cfg.foi;

%cfg.toi=toy;
cfg.toi       = -1:0.01:1;
%cfg.output         = 'pow';
cfg.output         = 'fourier';

cfg.keeptrials = 'yes';


freq = ft_freqanalysis(cfg, data1);




%%
cfg = [];
%      cfg.channelcmb = {'PFC' 'Parietal'};
     cfg.method    = 'granger';
     granger    = ft_connectivityanalysis(cfg, freq);

end
