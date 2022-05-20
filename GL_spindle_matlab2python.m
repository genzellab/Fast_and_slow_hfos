% GL_spindle_matlab2python.

% Main code for exporting NREM epochs from Matlab to python to then run the YASA algorithm 
%for spindle detection.

%% Select condition and session folder.

clear variables
close all
dname=uigetdir([],'Select folder with Matlab data of Rats trial');
cd(dname)

%%
%Band pass filter design:
fn=1000; % New sampling frequency. 

%LPF 300 Hz:
Wn1=[0.3/(fn/2) 300/(fn/2)]; % Cutoff=320 Hz
[b2,a2] = butter(3,Wn1,'bandpass'); %Filter coefficients

%Brain region combinations.
yy={'PAR'};    
xx={'PFC'};  
%%
PAR=dir(strcat('*',yy{1},'*.mat'));
PAR=PAR.name;
PAR=load(PAR);
PAR=getfield(PAR,yy{1});
PAR=PAR.*(0.195);

%PFC
PFC=dir(strcat('*',xx{1},'*.mat'));
PFC=PFC.name;
PFC=load(PFC);
PFC=getfield(PFC,xx{1});
PFC=PFC.*(0.195);



A = dir('*states*.mat');
A={A.name};

cellfun(@load,A);


    %Convert signal to 1 sec epochs.
    e_t=1;
    e_samples=e_t*(1000); %fs=1kHz
    ch=length(PAR);
    nc=floor(ch/e_samples); %Number of epochs
    NC=[];
    NC2=[];
    
    for kk=1:nc    
      NC(:,kk)= PAR(1+e_samples*(kk-1):e_samples*kk);
      NC2(:,kk)= PFC(1+e_samples*(kk-1):e_samples*kk);
    end
    
    vec_bin=states;
    vec_bin(vec_bin~=3)=0;
    vec_bin(vec_bin==3)=1;
    v2=ConsecutiveOnes(vec_bin);
    
    v_index=find(v2~=0);
    v_values=v2(v2~=0);

for epoch_count=1:length(v_index)
v_par{epoch_count,1}=reshape(NC(:, v_index(epoch_count):v_index(epoch_count)+(v_values(1,epoch_count)-1)), [], 1);
v_pfc{epoch_count,1}=reshape(NC2(:, v_index(epoch_count):v_index(epoch_count)+(v_values(1,epoch_count)-1)), [], 1);
end 

%v_par and v_pfc: NREM epochs.
%Low pass filter frequencies below 300Hz
V_par=cellfun(@(equis) filtfilt(b2,a2,equis), v_par ,'UniformOutput',false);
V_pfc=cellfun(@(equis) filtfilt(b2,a2,equis), v_pfc ,'UniformOutput',false);


%% Save data for YASA inputs.

save(['YASA'  '_' yy{1} '.mat'],'V_par');
save(['YASA'  '_' xx{1} '.mat'],'V_pfc');

msgbox('Run GL_yasa_spindles.py on Python');

