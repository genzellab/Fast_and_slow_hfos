function [values_spec,n,RandStruct]=getval_spectra(P,Q,labelconditions2,label1,s,w,win_size,same_nr_types,N,random_hfo,rand_first_run,tr)

if same_nr_types==1
   n=N;    
else
    n=min([length(P.(labelconditions2{1}).(label1{w}){s}) length(P.(labelconditions2{2}).(label1{w}){s})...
        length(P.(labelconditions2{3}).(label1{w}){s}) length(P.(labelconditions2{4}).(label1{w}){s})]);
    
end

type_hfo{1}='cohfos';
type_hfo{2}='single';


for condition=1:length(labelconditions2)

    %Order ripples
    p=P.(labelconditions2{condition}).(label1{w}){s}; 
    q=Q.(labelconditions2{condition}).(label1{w}){s}; 
    % R=(cellfun(@(equis1) max(abs(hilbert(equis1(1,:)))),q));
    %R=(cellfun(@(equis1) max(abs(hilbert(equis1(1,121-50:121+50)))),q));
    if random_hfo==0
    R=(cellfun(@(equis1) max(abs(hilbert(equis1(1,151-25:151+25)))),q));
    [~,r]=sort(R,'descend');
    p=p(r);
    q=q(r);
    p=p(1:n);
    q=q(1:n);
    else
        
            if rand_first_run==1 %First time generating rand numbers and saving them
                 rand_vec=randperm(length(p),n);
                 RandStruct.([label1{w} '_' type_hfo{s}]).(labelconditions2{condition})=rand_vec;
            else %Load old rand numbers for repeatability 
                
                if condition==1
                     if same_nr_types==1
                         load(['rand_' label1{w} '_' type_hfo{s} '_' num2str(tr(2)) '_SN.mat'])
                     else
                         load(['rand_' label1{w} '_' type_hfo{s} '_' num2str(tr(2)) '.mat'])
                     end
                end
                
                rand_vec=RandStruct.([label1{w} '_' type_hfo{s}]).(labelconditions2{condition});
                 
            end
            
        p=p(rand_vec);
        q=q(rand_vec);
    end


    %Max 1000 ripples.
    if length(q)>1000
        q=q(1:1000);
        p=p(1:1000);
    end

    if w==3 %PAR-centered ripples.
         p=cellfun(@(equis1) flip(equis1),p,'UniformOutput',false);
         q=cellfun(@(equis1) flip(equis1),q,'UniformOutput',false);
    end
    
    ro=150;

    toy=[-.1:.001:.1];
    freq=time_frequency(q,create_timecell(ro,length(q)),[100:1:300],[],toy);

    for j=1:3
    
    %Compute mean power value of window for different frequency bands
    [ndam,ndam2,ndam3,ndam4]=small_window(freq,j,win_size);
    
    
    Ndam(j)=ndam;
    Ndam2(j)=ndam2;
    Ndam3(j)=ndam3;
    Ndam4(j)=ndam4;

    end

    values_spec.(labelconditions2{condition})=[Ndam;Ndam2;Ndam3;Ndam4];
end

if random_hfo==1
     if same_nr_types==1
         save(['rand_' label1{w} '_' type_hfo{s} '_' num2str(tr(2)) '_SN.mat'],'RandStruct')
     else
         save(['rand_' label1{w} '_' type_hfo{s} '_' num2str(tr(2))  '.mat'],'RandStruct')
     end
end
                 
end