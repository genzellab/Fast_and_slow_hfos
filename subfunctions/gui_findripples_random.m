function [ripple,RipFreq,rip_duration,Mx,timeasleep,sig,Ex,Sx,ripple_multiplets,RipFreq_multiplets,rip_duration_multiplets,sig_multiplets,M_multiplets,Mr]=gui_findripples_random(CORTEX,states,xx,tr,multiplets,fn)
    %Band pass filter design:
    Wn1=[100/(fn/2) 300/(fn/2)]; % Cutoff=100-300 Hz
    [b1,a1] = butter(3,Wn1,'bandpass'); %Filter coefficients

    Wn1=[320/(fn/2)]; % Cutoff=320 Hz
    [b2,a2] = butter(3,Wn1); %Filter coefficients
%Convert signal to 1 sec epochs.
        e_t=1;
        e_samples=e_t*(fn); %fs=1kHz
        ch=length(CORTEX);
        nc=floor(ch/e_samples); %Number of epochsw
        NC=[];
        for kk=1:nc
          NC(:,kk)= CORTEX(1+e_samples*(kk-1):e_samples*kk);
        end
        vec_bin=states;
        vec_bin(vec_bin~=3)=0;
        vec_bin(vec_bin==3)=1;
        %Cluster one values:
        v2=ConsecutiveOnes(vec_bin);
        v_index=find(v2~=0);
        v_values=v2(v2~=0);
    for epoch_count=1:length(v_index)
    v{epoch_count,1}=reshape(NC(:, v_index(epoch_count):v_index(epoch_count)+(v_values(1,epoch_count)-1)), [], 1);
    end
    V=cellfun(@(equis) filtfilt(b2,a2,equis), v ,'UniformOutput',false);
    Mono=cellfun(@(equis) filtfilt(b1,a1,equis), V ,'UniformOutput',false);
    %Total amount of NREM time:
    timeasleep=sum(cellfun('length',V))*(1/fn)/60; % In minutes
    signal2=cellfun(@(equis) times((1/0.195), equis)  ,Mono,'UniformOutput',false);
    % ti=cellfun(@(equis) linspace(0, length(equis)-1,length(equis))*(1/fn) ,signal2,'UniformOutput',false);
    ti=cellfun(@(equis) reshape(linspace(0, length(equis)-1,length(equis))*(1/fn),[],1) ,signal2,'UniformOutput',false);

    %Find ripples or HFOs
    if strcmp(xx{1},'HPC')
    [Sx,Ex,Mx] =cellfun(@(equis1,equis2) findRipples(equis1, equis2, tr(1), (tr(1))*(1/2), fn ), signal2,ti,'UniformOutput',false);
    else
    [Sx,Ex,Mx] =cellfun(@(equis1,equis2) findHFOs(equis1, equis2, tr(2), (tr(2))*(1/2), [] ), signal2,ti,'UniformOutput',false);
    end
    
    %Generate 1000 shuffled timestamps of the events detected.
    for r=1:1000
        ti_rand=cellfun(@(equis) equis(randperm(size(equis, 1))),ti,'UniformOutput',false);
        Mr.(['Field_' num2str(r)])=cellfun(@(equis1,equis2,equis3) equis3(find(sum(equis1==equis2,2))).', ti,Mx,ti_rand,'UniformOutput',false );
        r
    end
    
    
    %Multiplets detection
    for l=1:length(Mx)
         hfo_sequence=ConsecutiveOnes(diff(Mx{l})<=0.300);

         for ll=1:length(multiplets)
             eval([multiplets{ll} '=(hfo_sequence=='  num2str(ll-1) ');'])
             cont=1;
             M_multiplets.(multiplets{ll}){l}=[];
             while cont<=ll
                 eval(['Mx_' multiplets{ll} '_' num2str(cont) '{l}=Mx{l}(find(' multiplets{ll} ')+(cont-1));'])
                 Mx_multiplets.(multiplets{ll}).(strcat('m_',num2str(cont))){l}=eval(['Mx_' multiplets{ll} '_' num2str(cont) '{l}']);
                 M_multiplets.(multiplets{ll}){l}=eval(['sort([M_multiplets.(multiplets{ll}){l} ' ' Mx_' multiplets{ll} '_' num2str(cont) '{l}])']); % Combined consecutive multiplets    
                 cont=cont+1;
             end
         end
    end
    
        %Multiplets detection
    for l=1:length(Mx)
         hfo_sequence=ConsecutiveOnes(diff(Mx{l})<=0.300);

         for ll=1:length(multiplets)
             eval([multiplets{ll} '=(hfo_sequence=='  num2str(ll-1) ');'])
             eval(['Sx_' multiplets{ll} '_1{l}=Sx{l}(find(' multiplets{ll} '));'])
             eval(['Ex_' multiplets{ll} '_1{l}=Ex{l}(find(' multiplets{ll} '));'])

         end
    end
    
    for l=1:length(Sx)
         sig{l}=getsignal(Sx,Ex,ti,Mono,l);
        for ll=1:length(multiplets)
            eval(['sig_' multiplets{ll} '_1{l}=getsignal(Sx_' multiplets{ll} '_1,Ex_' multiplets{ll} '_1,ti,Mono,l);'])
                        eval(strcat('sig_',multiplets{ll},'_1=sig_',multiplets{ll},'_1.'';'))
              eval(['sig_multiplets.(' 'multiplets{ll}' ')=sig_' multiplets{ll},'_1.'';'])  
              eval(['clear' ' ' 'sig_' multiplets{ll} '_1'])
        end
    end

    sig=sig.';

    [ripple, RipFreq,rip_duration]=hfo_count_freq_duration(Sx,Ex,timeasleep);
    
    for ll=1:length(multiplets)
        eval(['[ripple_multiplets.(' 'multiplets{ll}' '), RipFreq_multiplets.(' 'multiplets{ll}' '),rip_duration_multiplets.(' 'multiplets{ll}' ')]=hfo_count_freq_duration(Sx_' multiplets{ll} '_1,Ex_' multiplets{ll} '_1,timeasleep);']);
    end

end