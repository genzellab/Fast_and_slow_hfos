function [IND3]=gc_freqbands(gran,conditional,method)
% conditional=1 if conditional granger is to be used.

input_vec=gran.freq; %Get frequencies.
if strcmp(method,'granger') 
    input_gc=gran.grangerspctrm;
end
if strcmp(method,'psi')
    input_gc=gran.psispctrm;
    % input_gc(input_gc<0)=0; %ignore negative values.
end

%Frequency bands in which to inspect granger.
gc_bands=struct;
gc_bands.so=[0.01 4];
%gc_bands.delta=[2 4];
gc_bands.theta=[4 8];
gc_bands.spindle_beta=[10 20];
gc_bands.swr=[100 250];
% gc_bands.gamma_low=[30 50];
% gc_bands.gamma_high=[50 100];
% gc_bands.swr_slow=[100 150];
% gc_bands.swr_middle=[150 200];
% gc_bands.swr_fast=[200 250];
%Newest
% gc_bands.whole=[0 300]; %Whole spectrum
gc_bands.ZeroTwenty=[0 20];
gc_bands.TwentyThreehundred=[20 300];



GC_bands=struct2cell(gc_bands);
%Initialize variables
C=[];
c=[];
IND=[];
IND2=[];

    for j=1:length(fieldnames(gc_bands))
        freqlim=GC_bands{j};

        %Find index belonging to freq range.
        i1=max(find(input_vec<=freqlim(1)));
        i2=min(find(input_vec>=freqlim(2)));
        c{j}=input_vec(1,i1:i2); %Frequencies
        C(j,:)=[i1 i2]; %Index

        if conditional==0 %Not conditional
            IND{j}=input_gc(:,:,i1:i2);
            %Mean granger for the corresponding frequency band.
            IND2{j}=mean(input_gc(:,:,i1:i2),3); 
            %Different brain areas combinations.
            IND3(:,j)=[IND2{j}(1,2); IND2{j}(2,1); IND2{j}(1,3); IND2{j}(3,1); IND2{j}(2,3); IND2{j}(3,2)];
        else
            IND{j}=input_gc(:,i1:i2);
            %Mean granger for the corresponding frequency band.
            IND2{j}=mean(input_gc(:,i1:i2),2);
            IND3(:,:,j)=[IND2{j}(2) IND2{j}(1); IND2{j}(4) IND2{j}(3); IND2{j}(6) IND2{j}(5)];        
        end
    end
end