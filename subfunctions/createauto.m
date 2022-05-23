function [granger]=createauto(data1,ord,condition)
%Parametric model Order 20
cfg         = [];
cfg.order   = ord; %Model order
cfg.toolbox = 'bsmart';
mdata       = ft_mvaranalysis(cfg, data1);

cfg        = [];
cfg.method = 'mvar';
mfreq      = ft_freqanalysis(cfg, mdata);

cfg           = [];
cfg.method    = 'granger';

if strcmp(condition,'yes') %Conditional?
    cfg.granger.feedback    = 'yes';
    cfg.granger    = [];
    cfg.granger.conditional    = 'yes';   
end

granger       = ft_connectivityanalysis(cfg, mfreq);

end