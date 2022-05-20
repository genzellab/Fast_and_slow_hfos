function [Mono,V]=swr_preprocessing(CORTEX,states,fn)
%Band pass filter design:
    Wn1=[100/(fn/2) 300/(fn/2)]; % Cutoff=100-300 Hz
    [b1,a1] = butter(3,Wn1,'bandpass'); %Filter coefficients
    Wn1=[320/(fn/2)]; % Cutoff=320 Hz
    [b2,a2] = butter(3,Wn1); %Filter coefficients
%Convert signal to 1 sec epochs.
        e_t=1;
        e_samples=e_t*(fn); %fs=1kHz
        ch=length(CORTEX);
        nc=floor(ch/e_samples); %Number of epochs
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
    V=cellfun(@(equis) filtfilt(b2,a2,equis), v ,'UniformOutput',false); %Wideband
    Mono=cellfun(@(equis) filtfilt(b1,a1,equis), V ,'UniformOutput',false); %Bandpass
end