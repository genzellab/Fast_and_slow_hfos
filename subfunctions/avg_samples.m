function[P2]=avg_samples(p,timecell)
ft_data1 = [];
ft_data1.fsample = 1000;
ft_data1.trial = p; 

ft_data1.time = (timecell);

ft_data1.label = {'HPC'; 'PAR'; 'PFC'};

cfg = [];
    
   avg= ft_timelockanalysis(cfg,ft_data1);
   P2=avg.avg;

end