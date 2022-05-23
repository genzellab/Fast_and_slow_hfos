function [v1]=co_hfo_get_sample(cohf_mx_hpc,Cohfos1)

v1=[];
    for j=1:length(Cohfos1)
        
        n=find(cohf_mx_hpc==Cohfos1(j));
        v1=[v1 n];

    end

end