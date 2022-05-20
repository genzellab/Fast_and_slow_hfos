function [granger]=createauto_cond_multivariate(data1,ord)
%Parametric model Order 20
%Multivariate option set as on. (testing purposes)
cfg         = [];
cfg.order   = ord;
cfg.toolbox = 'bsmart';
mdata       = ft_mvaranalysis(cfg, data1);

cfg        = [];
cfg.method = 'mvar';
mfreq      = ft_freqanalysis(cfg, mdata);

cfg           = [];
cfg.method    = 'granger';

    cfg.granger.conditional    = 'yes';   
    
cfg.channelcmb  = {data1.label{1}, data1.label{2}, data1.label{3}};
cfg.granger.sfmethod = 'multivariate';

cfg.granger.block(1).name   = mfreq.label{1};
cfg.granger.block(1).label  = mfreq.label(1);
cfg.granger.block(2).name   = mfreq.label{2};
cfg.granger.block(2).label  = mfreq.label(2);
cfg.granger.block(3).name   = mfreq.label{3};
cfg.granger.block(3).label  = mfreq.label(3);


granger       = ft_connectivityanalysis(cfg, mfreq);

end