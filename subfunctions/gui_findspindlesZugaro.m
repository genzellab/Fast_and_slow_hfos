function [spindle,SpinFreq,spin_duration,Mx,timeasleep,sig,Ex,Sx]=gui_findspindlesZugaro(CORTEX,states,xx,multiplets,fn)
    %Band pass filter design:
    Wn1=[0.3/(fn/2) 300/(fn/2)]; 
    [b2,a2] = butter(3,Wn1); %0.3 to 300Hz
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
        %Convert to 1 if NREM.
        vec_bin(vec_bin~=3)=0;
        vec_bin(vec_bin==3)=1;
        %Cluster one values:
        v2=ConsecutiveOnes(vec_bin);
        v_index=find(v2~=0);
        v_values=v2(v2~=0);
    for epoch_count=1:length(v_index)
    v{epoch_count,1}=reshape(NC(:, v_index(epoch_count):v_index(epoch_count)+(v_values(1,epoch_count)-1)), [], 1);
    end
    V=cellfun(@(equis) filtfilt(b2,a2,equis), v ,'UniformOutput',false); %0.3 to 300Hz
    

    Wn1=[9/(fn/2) 20/(fn/2)]; % Cutoff=9-20 Hz
    [b1,a1] = butter(3,Wn1,'bandpass'); %Filter coefficients
    Mono=cellfun(@(equis) filtfilt(b1,a1,equis), V ,'UniformOutput',false); %Regular 9-20Hz bandpassed for sig variable.
    
    Concat_mono=vertcat(Mono{:});
      
    %Total amount of NREM time:
    timeasleep=sum(cellfun('length',V))*(1/fn)/60; % In minutes
    ti=cellfun(@(equis) reshape(linspace(0, length(equis)-1,length(equis))*(1/fn),[],1) ,Mono,'UniformOutput',false);
    
    %% Finding number of spindles in the dataset
    ti_cont=(1:length(Concat_mono))./1000;    
    Concat_input = [ti_cont' Concat_mono];
    spindles=FindSpindlesNayanika(Concat_input, 'durations', [0.5 2000], 'peak', 4);
    
    %finding number of spindles per epoch
    
    %duration of each epoch of non-rem sleep
    duration_epoch=cellfun(@length, ti)/1000; 
    duration_epoch_cumsum=cumsum(duration_epoch);
    
    for f=1:length(duration_epoch_cumsum)
        
        if(f==1)
            vec=find(spindles(:,1)>=0 & spindles(:,3)<=duration_epoch_cumsum(f));
        elseif(f>1)
            vec=find(spindles(:,1)>=duration_epoch_cumsum(f-1) & spindles(:,3)<=duration_epoch_cumsum(f));
        end
        
        spindle_per_epoch(f)=length(vec);
        % assigning each spindle an epoch number
        spindles(vec,5)=f;
        %assigning each spindle a local start and stop time wrt the epoch
        if(f==1)
            spindles(vec, 6:8)=spindles(vec,1:3);
        elseif(f>1)
            spindles(vec, 6:8)=spindles(vec,1:3)-duration_epoch_cumsum(f-1);
        end
         
        Sx(f,1)={spindles(vec,6)};
        Ex(f,1)={spindles(vec,8)};
        Mx(f,1)= {spindles(vec,7)};
    end
       
    %%
    % %% Find largest epoch.
    max_length=cellfun(@length,v);
    nrem_epoch=find(max_length==max(max_length)==1);
    nrem_epoch=nrem_epoch(1);
    
    % append all spindles together
    for l=1:length(Sx)
        sig{l}=getsignal(Sx,Ex,ti,Mono,l);
    end
    
    % Find # of spindles, spindle frequency, median spindle duration
    [spindle, SpinFreq,spin_duration]=hfo_count_freq_duration(Sx,Ex,timeasleep);
    
    SpinFreq=SpinFreq*60; %Spindles per minute.
end
