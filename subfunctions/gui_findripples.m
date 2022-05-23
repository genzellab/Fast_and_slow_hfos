function [ripple,RipFreq,rip_duration,Mx,timeasleep,sig,Ex,Sx,ripple_multiplets,RipFreq_multiplets,rip_duration_multiplets,sig_multiplets,M_multiplets,V,Mono,Sx_ti,Ex_ti,Mx_ti,v_rate]=gui_findripples(CORTEX,states,xx,tr,multiplets,fn)
    %Band pass filter design:
    Wn1=[100/(fn/2) 300/(fn/2)]; % Cutoff=100-300 Hz
    [b1,a1] = butter(3,Wn1,'bandpass'); %Filter coefficients
    Wn1=[320/(fn/2)]; % Cutoff=320 Hz
    [b2,a2] = butter(3,Wn1); %Filter coefficients
    ti_main=[0:length(CORTEX)-1]*(1/fn);
%Rearrange signal in 1 sec epochs.
        e_t=1;
        e_samples=e_t*(fn); %fs=1kHz
        ch=length(CORTEX);
        nc=floor(ch/e_samples); %Number of epochsw
        NC=[];
        NC_ti_main=[];
        for kk=1:nc
          NC(:,kk)= CORTEX(1+e_samples*(kk-1):e_samples*kk);
          NC_ti_main(:,kk)= ti_main(1+e_samples*(kk-1):e_samples*kk);
        end
        %Find if epoch is NREM (state=3)
        vec_bin=states;
        vec_bin(vec_bin~=3)=0;
        vec_bin(vec_bin==3)=1;
        %Cluster one values:
        v2=ConsecutiveOnes(vec_bin);
        v_index=find(v2~=0);
        v_values=v2(v2~=0);
    %Extract NREM epochs    
    for epoch_count=1:length(v_index)
    v{epoch_count,1}=reshape(NC(:, v_index(epoch_count):v_index(epoch_count)+(v_values(1,epoch_count)-1)), [], 1);
    v_ti{epoch_count,1}=reshape(NC_ti_main(:, v_index(epoch_count):v_index(epoch_count)+(v_values(1,epoch_count)-1)), [], 1);
    end
    
     % Splitting the 4 hour recording duration into 0.5 hour segments
     dur =0.5;
     labels=(0:1:length(states)-1);
     tet= [0:1800:((max(labels)-rem(max(labels),1800))+1800)]/60/60;
     if ceil(max(labels)/(30*60))>10
         xo % to check if the recording duration exceeds 5 hours
     end

     v_ti_sub{1,1} =cellfun(@(x) x(x/60/60<=tet(2)) , v_ti,'UniformOutput',false);
     ct=2;

     for ll = 1:(length(tet)-3)
         v_ti_sub{ct,1} =cellfun(@(x) x(x/60/60>(tet(ct)) & x/60/60<=(tet(ct)+dur)) , v_ti,'UniformOutput',false);
         ct=ct+1;
     end
     v_ti_sub{ct,1} =cellfun(@(x) x(x/60/60>tet(ct)) , v_ti,'UniformOutput',false);
    % Computing NREM duration in each hour
    for i=1:(length(tet)-1)
        v_rate(i,1)= sum(cell2mat(cellfun(@(x) (size(x,1)),v_ti_sub{i,1},'UniformOutput',false)))/(60*fn);
    end
    
    %Filter
    V=cellfun(@(equis) filtfilt(b2,a2,equis), v ,'UniformOutput',false);
    Mono=cellfun(@(equis) filtfilt(b1,a1,equis), V ,'UniformOutput',false);
    %Total amount of NREM time:
    timeasleep=sum(cellfun('length',V))*(1/fn)/60; % In minutes
    %NREM epochs without OpenEphys BitVolt factor
    signal2=cellfun(@(equis) times((1/0.195), equis)  ,Mono,'UniformOutput',false);

    %Timestamps for ripple/hfo detector
    ti=cellfun(@(equis) reshape(linspace(0, length(equis)-1,length(equis))*(1/fn),[],1) ,signal2,'UniformOutput',false);
    
    %If HPC find ripples, else find hfos.
    %Finds start, end and peak timestamp.
    if strcmp(xx{1},'HPC')
    [Sx,Ex,Mx] =cellfun(@(equis1,equis2) findRipples(equis1, equis2, tr(1), (tr(1))*(1/2), fn ), signal2,ti,'UniformOutput',false);
    else
    [Sx,Ex,Mx] =cellfun(@(equis1,equis2) findHFOs(equis1, equis2, tr(2), (tr(2))*(1/2), fn ), signal2,ti,'UniformOutput',false);
    end
    
    
    Sx_ti=cellfun(@(x,y,z) x(find(ismember(y*fn,z *fn))), v_ti, ti, Sx, 'UniformOutput', false );
    Ex_ti=cellfun(@(x,y,z) x(find(ismember(y*fn,z *fn))), v_ti, ti, Ex, 'UniformOutput', false );
    Mx_ti=cellfun(@(x,y,z) x(find(ismember(y*fn,z *fn))), v_ti, ti, Mx, 'UniformOutput', false );
    
    
    
%     % TIMESTAMPS WITH RESPECT TO RECORDING TIME.
%     %If HPC find ripples, else find hfos.
%     %Finds start, end and peak timestamp.
%     if strcmp(xx{1},'HPC')
%     [Sx_ti,Ex_ti,Mx_ti] =cellfun(@(equis1,equis2) findRipples(equis1, equis2, tr(1), (tr(1))*(1/2), fn ), signal2,v_ti,'UniformOutput',false);
%     else
%     [Sx_ti,Ex_ti,Mx_ti] =cellfun(@(equis1,equis2) findHFOs(equis1, equis2, tr(2), (tr(2))*(1/2), fn ), signal2,v_ti,'UniformOutput',false);
%     end

    
    %Multiplets detection
    for l=1:length(Mx)
         hfo_sequence=ConsecutiveOnes(diff(Mx{l})<=0.300);

         for ll=1:length(multiplets)
             eval([multiplets{ll} '=(hfo_sequence=='  num2str(ll-1) ');'])
             cont=1;
             M_multiplets.(multiplets{ll}){l}=[];
             while cont<=ll
                 %eval(['Sx_' multiplets{ll} '_' num2str(cont) '{l}=Sx{l}(find(' multiplets{ll} ')+(cont-1));'])
                 eval(['Mx_' multiplets{ll} '_' num2str(cont) '{l}=Mx{l}(find(' multiplets{ll} ')+(cont-1));'])
                 %eval(['Ex_' multiplets{ll} '_' num2str(cont) '{l}=Ex{l}(find(' multiplets{ll} ')+(cont-1));'])
                 Mx_multiplets.(multiplets{ll}).(strcat('m_',num2str(cont))){l}=eval(['Mx_' multiplets{ll} '_' num2str(cont) '{l}']);
                 M_multiplets.(multiplets{ll}){l}=eval(['sort([M_multiplets.(multiplets{ll}){l} ' ' Mx_' multiplets{ll} '_' num2str(cont) '{l}])']); % Combined consecutive multiplets    
%                   eval([  'clear' ' ' 'Mx_' multiplets{ll} '_' num2str(cont)])
                 cont=cont+1;
             end
         end
    end
    
    for l=1:length(Mx)
         hfo_sequence=ConsecutiveOnes(diff(Mx{l})<=0.300);

         for ll=1:length(multiplets)
             eval([multiplets{ll} '=(hfo_sequence=='  num2str(ll-1) ');'])
             eval(['Sx_' multiplets{ll} '_1{l}=Sx{l}(find(' multiplets{ll} '));'])
             eval(['Ex_' multiplets{ll} '_1{l}=Ex{l}(find(' multiplets{ll} '));'])

         end
    end

%Get traces of events detected    
    for l=1:length(Sx)
         sig{l}=getsignal(Sx,Ex,ti,Mono,l);
%         sig{l}=getsignal(Sx,Ex,ti,V,l);
        for ll=1:length(multiplets)
            eval(['sig_' multiplets{ll} '_1{l}=getsignal(Sx_' multiplets{ll} '_1,Ex_' multiplets{ll} '_1,ti,Mono,l);'])
%             eval(['sig_multiplets.(' 'multiplets{ll}' ')=getsignal(Sx_' multiplets{ll} '_1,Ex_' multiplets{ll} '_1,ti,Mono,l);'])
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
