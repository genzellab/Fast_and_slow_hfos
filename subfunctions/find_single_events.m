function [Sig_cortex_single]=find_single_events(Mx_cortex,cohfos2, sig_cortex)

v2=cellfun(@(equis1,equis2) single_hfo_get_sample(equis1,equis2),Mx_cortex,cohfos2,'UniformOutput',false);

Sig_cortex_single=cellfun(@(equis1,equis2) equis1(equis2),sig_cortex,v2,'UniformOutput',false);
Sig_cortex_single=[Sig_cortex_single{:}];

% [single_mx_cortex_val,single_sx_cortex_val]=cellfun(@(equis1,equis2,equis3) single_hfos_mx(equis1,equis2,equis3),cohfos2,Mx_cortex,Sx_cortex,'UniformOutput',false);
% single_mx_cortex_val=[single_mx_cortex_val{:}];
% single_sx_cortex_val=[single_sx_cortex_val{:}];

end