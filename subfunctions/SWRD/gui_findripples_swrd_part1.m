function [si,Sx,Ex,Mx,ti,Mono,V,timeasleep]=gui_findripples_swrd_part1(CORTEX,states,xx,tr,fn,wa2,Rat)
    %Band pass filter design:
    Wn1=[100/(fn/2) 300/(fn/2)]; % Cutoff=100-300 Hz
    [b1,a1] = butter(3,Wn1,'bandpass'); %Filter coefficients
    %LPF 300 Hz:
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
    if isempty(wa2)
        signal3=cellfun(@(equis) remove_stim_peaks(equis,Rat)  ,signal2,'UniformOutput',false); 
    else
        signal3=cellfun(@(equis1,equis2) equis1.*equis2  ,signal2,wa2,'UniformOutput',false); %Remove periods with HPC stimulation peaks
        signal3=cellfun(@(equis) remove_stim_peaks(equis,Rat)  ,signal3,'UniformOutput',false); %Remove large peaks in PPC.
    end
    
    % ti=cellfun(@(equis) linspace(0, length(equis)-1,length(equis))*(1/fn) ,signal2,'UniformOutput',false);
    ti=cellfun(@(equis) reshape(linspace(0, length(equis)-1,length(equis))*(1/fn),[],1) ,signal2,'UniformOutput',false);
    %xo
    if strcmp(xx{1},'HPC')
    [Sx,Ex,Mx] =cellfun(@(equis1,equis2) findRipples(equis1, equis2, tr(1), (tr(1))*(1/2), [] ), signal3,ti,'UniformOutput',false);
    else
    [Sx,Ex,Mx] =cellfun(@(equis1,equis2) findHFOs(equis1, equis2, tr(2), (tr(2))*(1/2), [] ), signal3,ti,'UniformOutput',false);
    end
   %%
    for l=1:length(Sx)
         sig{l}=getsignal(Sx,Ex,ti,Mono,l);
    end
   sig=sig.';
   
   %%
        
si=sig(~cellfun('isempty',sig));
si=[si{:}];
%%

end