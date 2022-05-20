function [data2,wa,wa2]=remove_stim_peaks(data,Rat)
    av=movmax(data, 100);
    if Rat==21 %|| Rat==24
    rng(100)
    else
    rng('Default')    
    end
    cluster = kmeans(av,2);   %input 1 must be column vector
    % Find indexes
    Idx = find(cluster == 2); 
    Idx2 = find(cluster == 1);
    meanVal=mean(data);
    %In case cluster numbers are inverted, make sure to use the one with
    %the largest absolute amplitude.
    
    if mean(av(Idx2(1))) > mean(av(Idx(1)))
        Idx=Idx2;
    end
    %%
    data2=data;
    %data2(Idx)=meanVal;     
    data2(Idx)=NaN;

   %Vector with ones and Nans to multiply with Cortical signal and omit periods of time with stimulation peaks. 
   wa=data2;
   wa(~isnan(wa))=1;  
   %%
veamos=wa;
veamos(isnan(veamos))=0;
cveamos=ConsecutiveOnes(veamos);
cveamos(cveamos<250)=0;
%%
wa2=nan(size(veamos));
f_cveamos=find(cveamos~=0);

    for cveamos_count=1:length(f_cveamos)
        wa2(f_cveamos(cveamos_count):f_cveamos(cveamos_count)+cveamos(f_cveamos(cveamos_count))-1)=1;
        %xo
    end

if sum(~isnan(wa2))==0 %When all elements are NaN
    wa2=zeros(size(wa2));
else
    wa2(isnan(wa2))=0;;
end

end