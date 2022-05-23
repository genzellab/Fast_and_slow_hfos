function [values_spec,n]=getval_spectra_All(P,Q,labelconditions2,label1,s,w,win_size)

n=([length(P.(labelconditions2{1}).(label1{w}){s}) length(P.(labelconditions2{2}).(label1{w}){s})...
        length(P.(labelconditions2{3}).(label1{w}){s}) length(P.(labelconditions2{4}).(label1{w}){s})]);

for condition=1:length(labelconditions2)

    %Order ripples
    p=P.(labelconditions2{condition}).(label1{w}){s}; 
    q=Q.(labelconditions2{condition}).(label1{w}){s}; 



    if w==3 %PAR-centered ripples.
         p=cellfun(@(equis1) flip(equis1),p,'UniformOutput',false);
         q=cellfun(@(equis1) flip(equis1),q,'UniformOutput',false);
    end
    
    ro=150;

    toy=[-.1:.001:.1];
    freq2=time_frequency(q,create_timecell(ro,length(q)),[100:1:300],[],toy);

    for j=1:3
    %Compute mean power value of window for different frequency bands    
    [ndam,ndam2,ndam3,ndam4]=small_window(freq2,j,win_size);
        
    Ndam(j)=ndam;
    Ndam2(j)=ndam2;
    Ndam3(j)=ndam3;
    Ndam4(j)=ndam4;

    end

    values_spec.(labelconditions2{condition})=[Ndam;Ndam2;Ndam3;Ndam4];
end

                 
end
