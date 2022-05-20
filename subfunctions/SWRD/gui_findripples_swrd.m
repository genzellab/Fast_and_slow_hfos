function [ripple2,RipFreq2,rip_duration,Mx,timeasleep,sig,Ex,Sx,ripple_multiplets,RipFreq_multiplets,rip_duration_multiplets,sig_multiplets,M_multiplets,V,Mono,wa,wa2]=gui_findripples_swrd(CORTEX,states,xx,tr,multiplets,fn,Rat)
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
    [signal3,wa,wa2]=cellfun(@(equis) remove_stim_peaks(equis,Rat)  ,signal2,'UniformOutput',false);
    
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
av=(cellfun(@(equis) max(abs(equis)),si));
index_false_positives=find(isoutlier(av)); %This index is with respect to si, which is the number of events, not epochs.
%%
% Correction on  Sx,Ex,Mx
Sx2=Sx;
counter=[];
Counter=[];
    for l=1:length(Sx2) %Epochs
        ind=Sx2{l};
        counter=[1: length(ind)]; %Amount of elements in given epoch
        if ~isempty(counter) 
            if l==1
                Counter=counter;
                c_count=counter;
            else
%                 c_count=max(Counter)+counter;
%                 Counter=[Counter max(Counter)+counter];
                if ~isempty(Counter)
                    c_count=max(Counter)+counter;
                    Counter=[Counter max(Counter)+counter];
                else
                    Counter=counter;
                    c_count=counter;    
                    
                end
            end
            detect_false=find(ismember(c_count,index_false_positives));
            if ~isempty(detect_false)
%                   xo
                stamps2remove=[findclosest([0: length(wa{l})-1]/fn,Sx{l}(detect_false)).' findclosest([0: length(wa{l})-1]/fn,Ex{l}(detect_false)).']; 
                for clu=1:size(stamps2remove,1)
                    copywa=wa{l};
                    copywa(stamps2remove(clu,1):stamps2remove(clu,2))=NaN; %Convert to NaN the false positives.
                    wa{l}=copywa;
                end
                Sx{l}(detect_false)=[];
                Mx{l}(detect_false)=[];
                Ex{l}(detect_false)=[];

            end
        end
%          sig{l}=getsignal(Sx,Ex,ti,Mono,l);
    end

clear sig
    for l=1:length(Sx)
         sig{l}=getsignal(Sx,Ex,ti,Mono,l);
    end

    sig=sig.';

if ~strcmp(xx{1},'HPC')
    
%% Second correction.

si=sig(~cellfun('isempty',sig));
si=[si{:}];
av=(cellfun(@(equis) max((equis)),si));
index_false_positives=find(isoutlier(av)); %This index is with respect to si, which is the number of events, not epochs.
% 
% % Correction on  Sx,Ex,Mx
Sx2=Sx;
counter=[];
Counter=[];
    for l=1:length(Sx2) %Epochs
        ind=Sx2{l};
        counter=[1: length(ind)]; %Amount of elements in given epoch
        if ~isempty(counter) 
            if l==1
                Counter=counter;
                c_count=counter;
            else
                c_count=max(Counter)+counter;
                Counter=[Counter max(Counter)+counter];
            end
            detect_false=find(ismember(c_count,index_false_positives));
            if ~isempty(detect_false)
               %    xo
                stamps2remove=[findclosest([0: length(wa{l})-1]/fn,Sx{l}(detect_false)).' findclosest([0: length(wa{l})-1]/fn,Ex{l}(detect_false)).']; 
                for clu=1:size(stamps2remove,1)
                    copywa=wa{l};
                    copywa(stamps2remove(clu,1):stamps2remove(clu,2))=NaN; %Convert to NaN the false positives.
                    wa{l}=copywa;
                end               
                Sx{l}(detect_false)=[];
                Mx{l}(detect_false)=[];
                Ex{l}(detect_false)=[];

            end
        end
%          sig{l}=getsignal(Sx,Ex,ti,Mono,l);
    end
% 
clear sig
    for l=1:length(Sx)
         sig{l}=getsignal(Sx,Ex,ti,Mono,l);
    end

    sig=sig.';
%% Third correction


si=sig(~cellfun('isempty',sig));
si=[si{:}];
av=(cellfun(@(equis) max(abs(equis)),si));
index_false_positives=find(isoutlier(av)); %This index is with respect to si, which is the number of events, not epochs.
% 
% % Correction on  Sx,Ex,Mx
Sx2=Sx;
counter=[];
Counter=[];
    for l=1:length(Sx2) %Epochs
        ind=Sx2{l};
        counter=[1: length(ind)]; %Amount of elements in given epoch
        if ~isempty(counter) 
            if l==1
                Counter=counter;
                c_count=counter;
            else
                c_count=max(Counter)+counter;
                Counter=[Counter max(Counter)+counter];
            end
            detect_false=find(ismember(c_count,index_false_positives));
            if ~isempty(detect_false)
               %    xo
                stamps2remove=[findclosest([0: length(wa{l})-1]/fn,Sx{l}(detect_false)).' findclosest([0: length(wa{l})-1]/fn,Ex{l}(detect_false)).']; 
                for clu=1:size(stamps2remove,1)
                    copywa=wa{l};
                    copywa(stamps2remove(clu,1):stamps2remove(clu,2))=NaN; %Convert to NaN the false positives.
                    wa{l}=copywa;
                end               
                Sx{l}(detect_false)=[];
                Mx{l}(detect_false)=[];
                Ex{l}(detect_false)=[];

            end
        end
%          sig{l}=getsignal(Sx,Ex,ti,Mono,l);
    end
% 
clear sig
    for l=1:length(Sx)
         sig{l}=getsignal(Sx,Ex,ti,Mono,l);
    end

    sig=sig.';
end
%% Multiplets
%     multiplets=[{'singlets'} {'doublets'} {'triplets'} {'quatruplets'} {'pentuplets'} {'sextuplets'} {'septuplets'} {'octuplets'} {'nonuplets'}];
    
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
    
        %Multiplets detection
    for l=1:length(Mx)
         hfo_sequence=ConsecutiveOnes(diff(Mx{l})<=0.300);

         for ll=1:length(multiplets)
             eval([multiplets{ll} '=(hfo_sequence=='  num2str(ll-1) ');'])
             eval(['Sx_' multiplets{ll} '_1{l}=Sx{l}(find(' multiplets{ll} '));'])
             eval(['Ex_' multiplets{ll} '_1{l}=Ex{l}(find(' multiplets{ll} '));'])

         end
    end
%     
%          
%          %Singlets
%          singlets=(hfo_sequence==0);
% %          Mx_singlets_1{l}=Mx{l}(find(singlets));
% %          Mx_singlets_2{l}=Mx{l}(find(singlets)+1);
%          Sx_singlets_1{l}=Sx{l}(find(singlets));
% %          Sx_singlets_2{l}=Sx{l}(find(singlets)+1);
%          Ex_singlets_1{l}=Ex{l}(find(singlets));
%          
%          %Douplets
%          douplets=(hfo_sequence==1);
% %          Mx_douplets_1{l}=Mx{l}(find(douplets));
% %          Mx_douplets_2{l}=Mx{l}(find(douplets)+1);
%          Sx_douplets_1{l}=Sx{l}(find(douplets));
% %          Sx_douplets_2{l}=Sx{l}(find(douplets)+1);
%          Ex_douplets_1{l}=Ex{l}(find(douplets));
% %          Ex_douplets_2{l}=Ex{l}(find(douplets)+1);
%          
%          %Triplets
%          triplets=(hfo_sequence==2);
% %          Mx_triplets_1{l}=Mx{l}(find(triplets));
% %          Mx_triplets_2{l}=Mx{l}(find(triplets)+1);
% %          Mx_triplets_3{l}=Mx{l}(find(triplets)+2);
%          Sx_triplets_1{l}=Sx{l}(find(triplets));
% %          Sx_triplets_2{l}=Sx{l}(find(triplets)+1);
% %          Sx_triplets_3{l}=Sx{l}(find(triplets)+2);
%          Ex_triplets_1{l}=Ex{l}(find(triplets));
% %          Ex_triplets_2{l}=Ex{l}(find(triplets)+1);
% %          Ex_triplets_3{l}=Ex{l}(find(triplets)+2);
%          
%          %Quadruplets
%          quadruplets=(hfo_sequence==3);
% %          Mx_quadruplets_1{l}=Mx{l}(find(quadruplets));
% %          Mx_quadruplets_2{l}=Mx{l}(find(quadruplets)+1);
% %          Mx_quadruplets_3{l}=Mx{l}(find(quadruplets)+2);
% %          Mx_quadruplets_4{l}=Mx{l}(find(quadruplets)+3);
%          
%          Sx_quadruplets_1{l}=Sx{l}(find(quadruplets));
% %          Sx_quadruplets_2{l}=Sx{l}(find(quadruplets)+1);
% %          Sx_quadruplets_3{l}=Sx{l}(find(quadruplets)+2);
% %          Sx_quadruplets_4{l}=Sx{l}(find(quadruplets)+3);
%          
%          Ex_quadruplets_1{l}=Ex{l}(find(quadruplets));
% %          Ex_quadruplets_2{l}=Ex{l}(find(quadruplets)+1);
% %          Ex_quadruplets_3{l}=Ex{l}(find(quadruplets)+2);
% %          Ex_quadruplets_4{l}=Ex{l}(find(quadruplets)+3);
%          
%          %Pentuplets
%          pentuplets=(hfo_sequence==4);
% %          Mx_pentuplets_1{l}=Mx{l}(find(pentuplets));
% %          Mx_pentuplets_2{l}=Mx{l}(find(pentuplets)+1);
% %          Mx_pentuplets_3{l}=Mx{l}(find(pentuplets)+2);
% %          Mx_pentuplets_4{l}=Mx{l}(find(pentuplets)+3);
% %          Mx_pentuplets_5{l}=Mx{l}(find(pentuplets)+4);
%          
%          
%          Sx_pentuplets_1{l}=Sx{l}(find(pentuplets));
% %          Sx_pentuplets_2{l}=Sx{l}(find(pentuplets)+1);
% %          Sx_pentuplets_3{l}=Sx{l}(find(pentuplets)+2);
% %          Sx_pentuplets_4{l}=Sx{l}(find(pentuplets)+3);
% %          Sx_pentuplets_5{l}=Sx{l}(find(pentuplets)+4);
% 
%          
%          Ex_pentuplets_1{l}=Ex{l}(find(pentuplets));
% %          Ex_pentuplets_2{l}=Ex{l}(find(pentuplets)+1);
% %          Ex_pentuplets_3{l}=Ex{l}(find(pentuplets)+2);
% %          Ex_pentuplets_4{l}=Ex{l}(find(pentuplets)+3);         
% %          Ex_pentuplets_5{l}=Ex{l}(find(pentuplets)+4);
%          
%          
%          %Sextuplets
%          sextuplets=(hfo_sequence==5);
% %          Mx_sextuplets_1{l}=Mx{l}(find(sextuplets));
% %          Mx_sextuplets_2{l}=Mx{l}(find(sextuplets)+1);
% %          Mx_sextuplets_3{l}=Mx{l}(find(sextuplets)+2);
% %          Mx_sextuplets_4{l}=Mx{l}(find(sextuplets)+3);
% %          Mx_sextuplets_5{l}=Mx{l}(find(sextuplets)+4);
% %          Mx_sextuplets_6{l}=Mx{l}(find(sextuplets)+5);
%                   
%          Sx_sextuplets_1{l}=Sx{l}(find(sextuplets));
% %          Sx_sextuplets_2{l}=Sx{l}(find(sextuplets)+1);
% %          Sx_sextuplets_3{l}=Sx{l}(find(sextuplets)+2);
% %          Sx_sextuplets_4{l}=Sx{l}(find(sextuplets)+3);
% %          Sx_sextuplets_5{l}=Sx{l}(find(sextuplets)+4);
% %          Sx_sextuplets_6{l}=Sx{l}(find(sextuplets)+5);
%          
%          Ex_sextuplets_1{l}=Ex{l}(find(sextuplets));
% %          Ex_sextuplets_2{l}=Ex{l}(find(sextuplets)+1);
% %          Ex_sextuplets_3{l}=Ex{l}(find(sextuplets)+2);
% %          Ex_sextuplets_4{l}=Ex{l}(find(sextuplets)+3);         
% %          Ex_sextuplets_5{l}=Ex{l}(find(sextuplets)+4);
% %          Ex_sextuplets_6{l}=Ex{l}(find(sextuplets)+5);  
% 
%          %Septuplets
%          septuplets=(hfo_sequence==6);
%          Sx_septuplets_1{l}=Sx{l}(find(septuplets));
%          Ex_septuplets_1{l}=Ex{l}(find(septuplets));
%          
%          %Octuplets
%          octuplets=(hfo_sequence==7);
%          Sx_octuplets_1{l}=Sx{l}(find(octuplets));
%          Ex_octuplets_1{l}=Ex{l}(find(octuplets));
% 
%          %Nonuplets
%          nonuplets=(hfo_sequence==8);
%          Sx_nonuplets_1{l}=Sx{l}(find(nonuplets));
%          Ex_nonuplets_1{l}=Ex{l}(find(nonuplets));
         
%     end    
    
    for l=1:length(Sx)
         sig{l}=getsignal(Sx,Ex,ti,Mono,l);
%        sig{l}=getsignal(Sx,Ex,ti,V,l);
        for ll=1:length(multiplets)
            eval(['sig_' multiplets{ll} '_1{l}=getsignal(Sx_' multiplets{ll} '_1,Ex_' multiplets{ll} '_1,ti,Mono,l);'])
%             eval(['sig_multiplets.(' 'multiplets{ll}' ')=getsignal(Sx_' multiplets{ll} '_1,Ex_' multiplets{ll} '_1,ti,Mono,l);'])
                        eval(strcat('sig_',multiplets{ll},'_1=sig_',multiplets{ll},'_1.'';'))
              eval(['sig_multiplets.(' 'multiplets{ll}' ')=sig_' multiplets{ll},'_1.'';'])  
              eval(['clear' ' ' 'sig_' multiplets{ll} '_1'])
        end
    end
    
    
    
    
%          sig_douplets_1{l}=getsignal(Sx_douplets_1,Ex_douplets_1,ti,Mono,l);
%          sig_triplets_1{l}=getsignal(Sx_triplets_1,Ex_triplets_1,ti,Mono,l);
%          sig_quadruplets_1{l}=getsignal(Sx_quadruplets_1,Ex_quadruplets_1,ti,Mono,l);
%          sig_pentuplets_1{l}=getsignal(Sx_pentuplets_1,Ex_pentuplets_1,ti,Mono,l);
%          sig_sextuplets_1{l}=getsignal(Sx_sextuplets_1,Ex_sextuplets_1,ti,Mono,l);         
%          sig_septuplets_1{l}=getsignal(Sx_septuplets_1,Ex_septuplets_1,ti,Mono,l);
%          sig_octuplets_1{l}=getsignal(Sx_octuplets_1,Ex_octuplets_1,ti,Mono,l); 
%          sig_nonuplets_1{l}=getsignal(Sx_nonuplets_1,Ex_nonuplets_1,ti,Mono,l);    

%     end



%%
    
%     sig_douplets_1=sig_douplets_1.';
%     sig_triplets_1=sig_triplets_1.';
%     sig_quadruplets_1=sig_quadruplets_1.';    
%     sig_pentuplets_1=sig_pentuplets_1.';    
%     sig_sextuplets_1=sig_sextuplets_1.';    
%     sig_septuplets_1=sig_septuplets_1.';    
%     sig_octuplets_1=sig_octuplets_1.';    
%     sig_nonuplets_1=sig_nonuplets_1.';    

    % [Sx,Ex,Mx] =cellfun(@(equis1,equis2) findHFOs(equis1, equis2, tr(2), (tr(2))*(1/2), [] ), signal2,ti,'UniformOutput',false);
%All HFOs
%     s=cellfun('length',Sx);
%     RipFreq2=sum(s)/(timeasleep*(60)); %RIpples per second.
%     ripple2=sum(s);
%     C = cellfun(@minus,Ex,Sx,'UniformOutput',false);
%     CC=([C{:}]);
%     rip_duration=median(CC);
    [ripple2, RipFreq2,rip_duration]=hfo_count_freq_duration(Sx,Ex,timeasleep);
    
    for ll=1:length(multiplets)
%        ripple_multiplets.(multiplets{ll})=1;
        eval(['[ripple_multiplets.(' 'multiplets{ll}' '), RipFreq_multiplets.(' 'multiplets{ll}' '),rip_duration_multiplets.(' 'multiplets{ll}' ')]=hfo_count_freq_duration(Sx_' multiplets{ll} '_1,Ex_' multiplets{ll} '_1,timeasleep);']);
    end

end