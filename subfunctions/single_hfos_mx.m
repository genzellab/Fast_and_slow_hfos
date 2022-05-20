
function [ach,ach2]=single_hfos_mx(cohfos1,ach,ach2)

    for k=1:length(cohfos1)
        
        ach2(find(ach==cohfos1(k)))=[];
        ach(find(ach==cohfos1(k)))=[];
   

    end

end