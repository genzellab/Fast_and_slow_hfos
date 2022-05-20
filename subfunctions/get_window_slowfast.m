function [p_cohfos_cortex,q_cohfos_cortex,p_single_cortex,q_single_cortex,use_me]=get_window_slowfast(Mx_cortex,Sx_cortex,Ex_cortex, cohfos2_g1,p_cortex,q_cortex,v2_g1)

%Cortical COHFOS
cohf_mx_cortex=Mx_cortex(~cellfun('isempty',cohfos2_g1));%Peak values cells where cortex cohfos were found.
cohf_sx_cortex=Sx_cortex(~cellfun('isempty',cohfos2_g1));%Start values cells where cortex cohfos were found.
cohf_ex_cortex=Ex_cortex(~cellfun('isempty',cohfos2_g1));%End values cells where cortex cohfos were found.

Cohfos2=cohfos2_g1(~cellfun('isempty',cohfos2_g1));

%Locate sample per cohfos
coh_samp_cortex= cellfun(@(equis1,equis2) co_hfo_get_sample(equis1,equis2),cohf_mx_cortex,Cohfos2,'UniformOutput',false);

%COHFOS windows
p_cohfos_cortex=p_cortex(~cellfun('isempty',cohfos2_g1));
p_cohfos_cortex=cellfun(@(equis1,equis2) equis1(equis2),p_cohfos_cortex,coh_samp_cortex,'UniformOutput',false);
p_cohfos_cortex=[p_cohfos_cortex{:}];
try
    use_me=~cellfun('isempty',p_cohfos_cortex);
catch
    warning('No ripples found.  Assigning an empty value.');
    use_me=[];
end
q_cohfos_cortex=q_cortex(~cellfun('isempty',cohfos2_g1));
q_cohfos_cortex=cellfun(@(equis1,equis2) equis1(equis2),q_cohfos_cortex,coh_samp_cortex,'UniformOutput',false);
q_cohfos_cortex=[q_cohfos_cortex{:}];
if ~isempty(p_cohfos_cortex)
p_cohfos_cortex=p_cohfos_cortex(~cellfun('isempty',p_cohfos_cortex));
q_cohfos_cortex=q_cohfos_cortex(~cellfun('isempty',q_cohfos_cortex));
else
p_cohfos_cortex=[];
q_cohfos_cortex=[];

end



%Single cortex windows
p_single_cortex=cellfun(@(equis1,equis2) equis1(equis2),p_cortex,v2_g1,'UniformOutput',false);
p_single_cortex=[p_single_cortex{:}];
q_single_cortex=cellfun(@(equis1,equis2) equis1(equis2),q_cortex,v2_g1,'UniformOutput',false);
q_single_cortex=[q_single_cortex{:}];
p_single_cortex=p_single_cortex(~cellfun('isempty',p_single_cortex));
q_single_cortex=q_single_cortex(~cellfun('isempty',q_single_cortex));

end