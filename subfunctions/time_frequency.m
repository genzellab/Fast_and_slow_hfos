function [freq]=time_frequency(q,timecell,freqrange,label,toy,pad_option)

ft_data1 = [];
ft_data1.fsample = 1000;
ft_data1.trial = q(1,1:end); 
ft_data1.time = (timecell(1,1:end));

ft_data1
ft_data1.label = {'HPC'; 'PAR'; 'PFC'};

% Compute spectrogram

cfg = [];
cfg.method = 'mtmconvol';
cfg.taper = 'dpss';
cfg.foi = freqrange;
cfg.t_ftimwin = .1 * ones(size(cfg.foi));
cfg.tapsmofrq = 10;
cfg.toi=toy;
cfg.keeptrials = 'yes';
cfg.output         = 'pow';

freq = ft_freqanalysis(cfg, ft_data1);

end