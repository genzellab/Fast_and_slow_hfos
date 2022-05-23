function [v2]=single_hfo_get_sample(Mx_hpc,cohfos1)

a=1:length(Mx_hpc);
v1=co_hfo_get_sample(Mx_hpc,cohfos1);
v2=a(~ismember(a,v1));

end