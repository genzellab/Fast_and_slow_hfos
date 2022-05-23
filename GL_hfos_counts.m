%% GL_hfos_counts

% Main code for detection of fast and slow, coocur and single hfos.
% Generates hfos counts from Figure 1 and 2 respectively.
% Requires 'load_me_first.mat' loaded first. /home/adrian/Documents/GitHub/CorticoHippocampal/Fast_and_slow_hfos
clear variables
cd('/home/adrian/Documents/GitHub/CorticoHippocampal/Fast_and_slow_hfos')
load('load_me_first.mat')
%% Find location of downsampled data

close all
% dname=uigetdir([],'Select folder with Matlab data containing all rats.');
%dname='/home/adrian/Dropbox/jukebox/Documents/Plusmaze_downsampled';
dname='/media/adrian/6aa1794c-0320-4096-a7df-00ab0ba946dc/Plusmaze_downsampled/Data_plusmaze';
cd(dname)

%cd('/home/adrian/Documents/Plusmaze_downsampled')
%/home/adrian/Dropbox/jukebox/Documents/Plusmaze_downsampled

%%
%Select rat ID
opts.Resize = 'on';
opts.WindowStyle = 'modal';
opts.Interpreter = 'tex';
prompt=strcat('\bf Select a rat#. Options:','{ }',num2str(rats));
answer = inputdlg(prompt,'Input',[2 30],{''},opts);
Rat=str2num(answer{1});
cd(num2str(Rat))
tr=getfield(T,strcat('Rat',num2str(Rat)));%Thresholds 
%%
xx={'PAR'}; %Posterior Parietal cortex used to detect hfos. 
fn=1000; %Sampling frequency

%% Get folder names
labelconditions=getfolder;
labelconditions=labelconditions.';
g=labelconditions;
%% Parameters 
multiplets=[{'singlets'} {'doublets'} {'triplets'} {'quatruplets'} {'pentuplets'} {'sextuplets'} {'septuplets'} {'octuplets'} {'nonuplets'}];
iii=1;
%%
f=waitbar(0,'Please wait...');
    for k=1:length(g)
        cd(g{k})
%xo
CORTEX=dir(strcat('*',xx{1},'*.mat'));
if isempty(CORTEX)
    g=g(~contains(g,g{k}));
    cd ..
    progress_bar(k,length(g),f)
    break
end
CORTEX=CORTEX.name;
CORTEX=load(CORTEX);
%CORTEX=CORTEX.CORTEX;
CORTEX=getfield(CORTEX,xx{1});
CORTEX=CORTEX.*(0.195); % OpenEphys BitVolt factor

%Load sleep scoring
A = dir('*states*.mat');
A={A.name};

if  ~isempty(A)
       cellfun(@load,A);
else
      error('No Scoring found')    
end
%xo
%Find hfos
[ripple,RipFreq,rip_duration,Mx_cortex,timeasleep,sig_cortex,Ex_cortex,Sx_cortex,...
  ripple_multiplets_cortex,RipFreq_multiplets_cortex,rip_duration_multiplets_cortex,sig_multiplets_cortex,~,~,~,Sx_cortex_ti,Ex_cortex_ti,Mx_cortex_ti,v_rate_CR{k} ...
  ]=gui_findripples(CORTEX,states,xx,tr,multiplets,fn);

[Mx_cortex_ti_1] =cellfun(@(x) x(x/60/60<=1) , Mx_cortex_ti,'UniformOutput',false);
[Mx_cortex_ti_2] =cellfun(@(x) x(x/60/60>1 & x/60/60<=2) , Mx_cortex_ti,'UniformOutput',false);
[Mx_cortex_ti_3] =cellfun(@(x) x(x/60/60>2 & x/60/60<=3) , Mx_cortex_ti,'UniformOutput',false);
[Mx_cortex_ti_4] =cellfun(@(x) x(x/60/60>3) , Mx_cortex_ti,'UniformOutput',false);

[Sx_cortex_ti_1] =cellfun(@(x) x(x/60/60<=1) , Sx_cortex_ti,'UniformOutput',false);
[Sx_cortex_ti_2] =cellfun(@(x) x(x/60/60>1 & x/60/60<=2) , Sx_cortex_ti,'UniformOutput',false);
[Sx_cortex_ti_3] =cellfun(@(x) x(x/60/60>2 & x/60/60<=3) , Sx_cortex_ti,'UniformOutput',false);
[Sx_cortex_ti_4] =cellfun(@(x) x(x/60/60>3) , Sx_cortex_ti,'UniformOutput',false);

%% Plot cortical ripples histogram across 4 hours.
figure()
xlabel('Time (Hours)','FontSize',12)
ylabel('Amount of ripples','FontSize',12)
title('Histogram of cortical ripples','FontSize',12)

vec_trans=(states==4);
vec_rem=(states==5);
vec_nrem=(states==3);
vec_wake=not(vec_trans) & not(vec_rem) & not(vec_nrem);
labels=(0:1:length(states)-1);

rip_times=cell2mat(Mx_cortex_ti);
hist_rip=histcounts(rip_times,[0:10: max(labels)+1]);
CR_30_count(k)={histcounts(rip_times,[0:1800:((max(labels)-rem(max(labels),1800))+1800)])};
if ceil(max(labels)/(30*60))>10
    xo % to check if the recording duration exceeds 5 hours
end

yu = 35; % ylim
hold on
stripes(vec_trans,0.2,labels/60/60,'b',yu)
stripes(vec_rem,0.2,labels/60/60,'r',yu)
stripes(vec_nrem,0.9,labels/60/60,'k',yu)


stem(linspace(0,max(labels)/60/60,length(hist_rip)),hist_rip,'filled','Color',[0.3010 0.7450 0.9330])
ylim([0 yu])
xlim([0 4])

cd ..
 printing(['Histograms_CorticalRipples_across_time' g{k}]);
close all
cd(g{k})

%%
%Get traces of events and store them
si=sig_cortex(~cellfun('isempty',sig_cortex));
si=[si{:}];
All_Par.( strrep(g{k},'-','_'))=si;

%Compute main features of hfos.
print_hist=0;
[x,y,z,~,~,~,l,p,si_mixed,th]=hfo_specs(si,timeasleep,print_hist,Rat,tr,fn);
if print_hist==1
    cd ..
    printing(['Histograms_Cortex_Count_' g{k}]);
    close all
    cd(g{k})
end

fi_cortex(k)=x; %Instantaneous frequency.
fa_cortex(k)=y; %Average frequency.
amp_cortex(k)=z; %Amplitude.
auc_cortex(k)=l; %Area under the curve.
p2p_cortex(k)=p; %Peak to trough distance
TH(k)=th;
%% Cortical HFOs
    hfos_cortex(k)=ripple;
    hfos_cortex_rate(k)=RipFreq;
    hfos_cortex_duration(k)=rip_duration;
    clear ripple RipFreq

%Multiplets    
for ll=1:3
   eval(['hfos_cortex_' multiplets{ll} '(k)=ripple_multiplets_cortex.' multiplets{ll} ';']) 
   eval(['hfos_cortex_rate_' multiplets{ll} '(k)=RipFreq_multiplets_cortex.' multiplets{ll} ';']) 
   eval(['hfos_cortex_duration_' multiplets{ll} '(k)=rip_duration_multiplets_cortex.' multiplets{ll} ';'])    
end
    
%% HPC     
HPC=dir(strcat('*','HPC','*.mat'));
HPC=HPC.name;
HPC=load(HPC);
HPC=getfield(HPC,'HPC');
HPC=HPC.*(0.195);
%xo
%Find ripples 
[ripple,RipFreq,rip_duration,Mx_hpc,timeasleep,sig_hpc,Ex_hpc,Sx_hpc,...
  ripple_multiplets_hpc,RipFreq_multiplets_hpc,rip_duration_multiplets_hpc,sig_multiplets_hpc,Mx_multiplets_hpc,~,~,Sx_ti,Ex_ti,Mx_ti,v_rate...    
  ]=gui_findripples(HPC,states,{'HPC'},tr,multiplets,fn);

[Mx_ti_1] =cellfun(@(x) x(x/60/60<=1) , Mx_ti,'UniformOutput',false);
[Mx_ti_2] =cellfun(@(x) x(x/60/60>1 & x/60/60<=2) , Mx_ti,'UniformOutput',false);
[Mx_ti_3] =cellfun(@(x) x(x/60/60>2 & x/60/60<=3) , Mx_ti,'UniformOutput',false);
[Mx_ti_4] =cellfun(@(x) x(x/60/60>3) , Mx_ti,'UniformOutput',false);



%% Plot hippocampal ripples histogram across 4 hours.
figure()
xlabel('Time (Hours)','FontSize',12)
ylabel('Amount of ripples','FontSize',12)
title('Histogram of hippocampal ripples','FontSize',12)

vec_trans=(states==4);
vec_rem=(states==5);
vec_nrem=(states==3);
vec_wake=not(vec_trans) & not(vec_rem) & not(vec_nrem);
labels=(0:1:length(states)-1);

rip_times=cell2mat(Mx_ti);
hist_rip=histcounts(rip_times,[0:10: max(labels)+1]);
HR_30_count(k)={histcounts(rip_times,[0:1800:((max(labels)-rem(max(labels),1800))+1800)])};
if ceil(max(labels)/(30*60))>10
    xo % to check if the recording duration exceeds 5 hours
end

yu = 35; % ylim

hold on
stripes(vec_trans,0.2,labels/60/60,'b',yu)
stripes(vec_rem,0.2,labels/60/60,'r',yu)
stripes(vec_nrem,0.9,labels/60/60,'k',yu)


stem(linspace(0,max(labels)/60/60,length(hist_rip)),hist_rip,'filled','Color',[0.3010 0.7450 0.9330])
ylim([0 yu])
xlim([0 4])

cd ..
 printing(['Histograms_HippocampalRipples_across_time' g{k}]);
close all
cd(g{k})


%%
%Ripple traces
si=sig_hpc(~cellfun('isempty',sig_hpc));
si=[si{:}];
All_HPC.( strrep(g{k},'-','_'))=si;

%Compute ripples main features
print_hist=0; %No printing
[x,y,z,~,~,~,l,p]=hfo_specs_hpc(si,timeasleep,print_hist,fn);
if print_hist==1
    cd ..
    printing(['Histograms_HPC_Probability_' g{k}]);
    close all
    cd(g{k})
end

fi_hpc(k)=x;
fa_hpc(k)=y;
amp_hpc(k)=z;
auc_hpc(k)=l;
p2p_hpc(k)=p;
% xo
%% HPC HFOs
hfos_hpc(k)=ripple;
hfos_hpc_rate(k)=RipFreq;
hfos_hpc_duration(k)=rip_duration;

hfos_hpc_perhour(k,1)=sum(cellfun('length',Mx_ti_1));
hfos_hpc_perhour(k,2)=sum(cellfun('length',Mx_ti_2));
hfos_hpc_perhour(k,3)=sum(cellfun('length',Mx_ti_3));
hfos_hpc_perhour(k,4)=sum(cellfun('length',Mx_ti_4));


%Multiplets    
for ll=1:length(multiplets)
   eval(['hfos_hpc_' multiplets{ll} '(k)=ripple_multiplets_hpc.' multiplets{ll} ';']) 
   eval(['hfos_hpc_rate_' multiplets{ll} '(k)=RipFreq_multiplets_hpc.' multiplets{ll} ';']) 
   eval(['hfos_hpc_duration_' multiplets{ll} '(k)=rip_duration_multiplets_hpc.' multiplets{ll} ';'])    
end
%xo
%% Coocurent hfos and ripples
[cohfos1,cohfos2]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),Mx_hpc,Mx_cortex,'UniformOutput',false);
%cohfos1: HPC timestamps.
%cohfos2: Parietal timestamps.
%Counts and rate:
cohfos_count(k)=sum(cellfun('length',cohfos1));
cohfos_rate(k)=sum(cellfun('length',cohfos1))/(timeasleep*(60));

%Cohfos per hour
[cohfos1_ti_1,cohfos2_ti_1]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),Mx_ti_1,Mx_cortex_ti_1,'UniformOutput',false);
[cohfos1_ti_2,cohfos2_ti_2]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),Mx_ti_2,Mx_cortex_ti_2,'UniformOutput',false);
[cohfos1_ti_3,cohfos2_ti_3]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),Mx_ti_3,Mx_cortex_ti_3,'UniformOutput',false);
[cohfos1_ti_4,cohfos2_ti_4]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),Mx_ti_4,Mx_cortex_ti_4,'UniformOutput',false);

cohfos_count_perhour(k,1)=sum(cellfun('length',cohfos1_ti_1));
cohfos_count_perhour(k,2)=sum(cellfun('length',cohfos1_ti_2));
cohfos_count_perhour(k,3)=sum(cellfun('length',cohfos1_ti_3));
cohfos_count_perhour(k,4)=sum(cellfun('length',cohfos1_ti_4));

cohfos_rate_perhour(k,1)=sum(cellfun('length',cohfos1_ti_1))/(v_rate(1)*(60));
cohfos_rate_perhour(k,2)=sum(cellfun('length',cohfos1_ti_2))/(v_rate(2)*(60));
cohfos_rate_perhour(k,3)=sum(cellfun('length',cohfos1_ti_3))/(v_rate(3)*(60));
cohfos_rate_perhour(k,4)=sum(cellfun('length',cohfos1_ti_4))/(v_rate(4)*(60));


%Multiplet cohfos
for ll=1:length(multiplets)
[cohfos1_multiplets.(multiplets{ll}),cohfos2_multiplets.(multiplets{ll})]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),Mx_multiplets_hpc.(multiplets{ll}).',Mx_cortex,'UniformOutput',false);
cohfos_count_multiplets.(multiplets{ll})(k)=sum(cellfun('length',cohfos1_multiplets.(multiplets{ll})));
cohfos_rate_multiplets.(multiplets{ll})(k)=sum(cellfun('length',cohfos1_multiplets.(multiplets{ll})))/(timeasleep*(60));
end

%% Mixed distribution (Average freq) coHFOs
%Find timestamps corresponding to g1 and g2 events.
%g1 and g2 stand for gaussian1 and gaussian2 (slow and fast hfos)
Mx_cortex_g1=Mx_cortex;
Mx_cortex_g2=Mx_cortex;



row=si_mixed.i1;
cont=0;
for ll=1:length(Mx_cortex)

    if ~isempty(Mx_cortex{ll})

        for lll=1:length(Mx_cortex{ll})
            cont=cont+1;
    %         xo

            if ~ismember(cont,row)
                Mx_cortex_g1{ll}(lll)=NaN;
            else
                Mx_cortex_g2{ll}(lll)=NaN;
            end

        end
         Mx_cortex_g1{ll}=Mx_cortex_g1{ll}(~isnan(Mx_cortex_g1{ll}));
         Mx_cortex_g2{ll}=Mx_cortex_g2{ll}(~isnan(Mx_cortex_g2{ll}));

    end

end
%Coocur slow and coocur fast
[cohfos1_g1,cohfos2_g1]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),Mx_hpc,Mx_cortex_g1,'UniformOutput',false);
[cohfos1_g2,cohfos2_g2]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),Mx_hpc,Mx_cortex_g2,'UniformOutput',false);

cohfos_count_g1(k)=sum(cellfun('length',cohfos1_g1));
cohfos_rate_g1(k)=sum(cellfun('length',cohfos1_g1))/(timeasleep*(60));

cohfos_count_g2(k)=sum(cellfun('length',cohfos1_g2));
cohfos_rate_g2(k)=sum(cellfun('length',cohfos1_g2))/(timeasleep*(60));

%% Find slow and fast per hour

Mx_cortex_g1_ti=Mx_cortex_ti;
Mx_cortex_g2_ti=Mx_cortex_ti;



row=si_mixed.i1;
cont=0;
for ll=1:length(Mx_cortex)

    if ~isempty(Mx_cortex{ll})

        for lll=1:length(Mx_cortex{ll})
            cont=cont+1;
    %         xo

            if ~ismember(cont,row)
                Mx_cortex_g1_ti{ll}(lll)=NaN;
            else
                Mx_cortex_g2_ti{ll}(lll)=NaN;
            end

        end
         Mx_cortex_g1_ti{ll}=Mx_cortex_g1_ti{ll}(~isnan(Mx_cortex_g1_ti{ll}));
         Mx_cortex_g2_ti{ll}=Mx_cortex_g2_ti{ll}(~isnan(Mx_cortex_g2_ti{ll}));

    end

end
%%
[Mx_cortex_g1_ti_1] =cellfun(@(x) x(x/60/60<=1) , Mx_cortex_g1_ti,'UniformOutput',false);
[Mx_cortex_g1_ti_2] =cellfun(@(x) x(x/60/60>1 & x/60/60<=2) , Mx_cortex_g1_ti,'UniformOutput',false);
[Mx_cortex_g1_ti_3] =cellfun(@(x) x(x/60/60>2 & x/60/60<=3) , Mx_cortex_g1_ti,'UniformOutput',false);
[Mx_cortex_g1_ti_4] =cellfun(@(x) x(x/60/60>3) , Mx_cortex_g1_ti,'UniformOutput',false);

[Mx_cortex_g2_ti_1] =cellfun(@(x) x(x/60/60<=1) , Mx_cortex_g2_ti,'UniformOutput',false);
[Mx_cortex_g2_ti_2] =cellfun(@(x) x(x/60/60>1 & x/60/60<=2) , Mx_cortex_g2_ti,'UniformOutput',false);
[Mx_cortex_g2_ti_3] =cellfun(@(x) x(x/60/60>2 & x/60/60<=3) , Mx_cortex_g2_ti,'UniformOutput',false);
[Mx_cortex_g2_ti_4] =cellfun(@(x) x(x/60/60>3) , Mx_cortex_g2_ti,'UniformOutput',false);



%Coocur slow and coocur fast
% HOUR 1
[cohfos1_g1_1,cohfos2_g1_1]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),Mx_ti_1,Mx_cortex_g1_ti_1,'UniformOutput',false);
[cohfos1_g2_1,cohfos2_g2_1]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),Mx_ti_1,Mx_cortex_g2_ti_1,'UniformOutput',false);

cohfos_count_g1_ti(k,1)=sum(cellfun('length',cohfos1_g1_1));
cohfos_rate_g1_ti(k,1)=sum(cellfun('length',cohfos1_g1_1))/(v_rate(1)*(60));

cohfos_count_g2_ti(k,1)=sum(cellfun('length',cohfos1_g2_1));
cohfos_rate_g2_ti(k,1)=sum(cellfun('length',cohfos1_g2_1))/(v_rate(1)*(60));

% HOUR 2
[cohfos1_g1_2,cohfos2_g1_2]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),Mx_ti_2,Mx_cortex_g1_ti_2,'UniformOutput',false);
[cohfos1_g2_2,cohfos2_g2_2]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),Mx_ti_2,Mx_cortex_g2_ti_2,'UniformOutput',false);

cohfos_count_g1_ti(k,2)=sum(cellfun('length',cohfos1_g1_2));
cohfos_rate_g1_ti(k,2)=sum(cellfun('length',cohfos1_g1_2))/(v_rate(2)*(60));

cohfos_count_g2_ti(k,2)=sum(cellfun('length',cohfos1_g2_2));
cohfos_rate_g2_ti(k,2)=sum(cellfun('length',cohfos1_g2_2))/(v_rate(2)*(60));

% HOUR 3
[cohfos1_g1_3,cohfos2_g1_3]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),Mx_ti_3,Mx_cortex_g1_ti_3,'UniformOutput',false);
[cohfos1_g2_3,cohfos2_g2_3]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),Mx_ti_3,Mx_cortex_g2_ti_3,'UniformOutput',false);

cohfos_count_g1_ti(k,3)=sum(cellfun('length',cohfos1_g1_3));
cohfos_rate_g1_ti(k,3)=sum(cellfun('length',cohfos1_g1_3))/(v_rate(3)*(60));

cohfos_count_g2_ti(k,3)=sum(cellfun('length',cohfos1_g2_3));
cohfos_rate_g2_ti(k,3)=sum(cellfun('length',cohfos1_g2_3))/(v_rate(3)*(60));

% HOUR 4
[cohfos1_g1_4,cohfos2_g1_4]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),Mx_ti_4,Mx_cortex_g1_ti_4,'UniformOutput',false);
[cohfos1_g2_4,cohfos2_g2_4]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),Mx_ti_4,Mx_cortex_g2_ti_4,'UniformOutput',false);

cohfos_count_g1_ti(k,4)=sum(cellfun('length',cohfos1_g1_4));
cohfos_rate_g1_ti(k,4)=sum(cellfun('length',cohfos1_g1_4))/(v_rate(4)*(60));

cohfos_count_g2_ti(k,4)=sum(cellfun('length',cohfos1_g2_4));
cohfos_rate_g2_ti(k,4)=sum(cellfun('length',cohfos1_g2_4))/(v_rate(4)*(60));
%xo
%-----------
%Find single fast and slow events
v2_g1=cellfun(@(equis1,equis2) single_hfo_get_sample(equis1,equis2),Mx_cortex_g1,cohfos2_g1,'UniformOutput',false);
singles_count_g1(k)=sum(cellfun('length',v2_g1));
singles_rate_g1(k)=sum(cellfun('length',v2_g1))/(timeasleep*(60));


v2_g2=cellfun(@(equis1,equis2) single_hfo_get_sample(equis1,equis2),Mx_cortex_g2,cohfos2_g2,'UniformOutput',false);
singles_count_g2(k)=sum(cellfun('length',v2_g2));
singles_rate_g2(k)=sum(cellfun('length',v2_g2))/(timeasleep*(60));
%------------
%% HOUR 1
%Find single fast and slow events
v2_g1_1=cellfun(@(equis1,equis2) single_hfo_get_sample(equis1,equis2),Mx_cortex_g1_ti_1,cohfos2_g1_1,'UniformOutput',false);
singles_count_g1_ti(k,1)=sum(cellfun('length',v2_g1_1));
singles_rate_g1_ti(k,1)=sum(cellfun('length',v2_g1_1))/(v_rate(1)*(60));


v2_g2_1=cellfun(@(equis1,equis2) single_hfo_get_sample(equis1,equis2),Mx_cortex_g2_ti_1,cohfos2_g2_1,'UniformOutput',false);
singles_count_g2_ti(k,1)=sum(cellfun('length',v2_g2_1));
singles_rate_g2_ti(k,1)=sum(cellfun('length',v2_g2_1))/(v_rate(1)*(60));

%% HOUR 2
%Find single fast and slow events
v2_g1_2=cellfun(@(equis1,equis2) single_hfo_get_sample(equis1,equis2),Mx_cortex_g1_ti_2,cohfos2_g1_2,'UniformOutput',false);
singles_count_g1_ti(k,2)=sum(cellfun('length',v2_g1_2));
singles_rate_g1_ti(k,2)=sum(cellfun('length',v2_g1_2))/(v_rate(2)*(60));


v2_g2_2=cellfun(@(equis1,equis2) single_hfo_get_sample(equis1,equis2),Mx_cortex_g2_ti_2,cohfos2_g2_2,'UniformOutput',false);
singles_count_g2_ti(k,2)=sum(cellfun('length',v2_g2_2));
singles_rate_g2_ti(k,2)=sum(cellfun('length',v2_g2_2))/(v_rate(2)*(60));

%% HOUR 3
%Find single fast and slow events
v2_g1_3=cellfun(@(equis1,equis2) single_hfo_get_sample(equis1,equis2),Mx_cortex_g1_ti_3,cohfos2_g1_3,'UniformOutput',false);
singles_count_g1_ti(k,3)=sum(cellfun('length',v2_g1_3));
singles_rate_g1_ti(k,3)=sum(cellfun('length',v2_g1_3))/(v_rate(3)*(60));


v2_g2_3=cellfun(@(equis1,equis2) single_hfo_get_sample(equis1,equis2),Mx_cortex_g2_ti_3,cohfos2_g2_3,'UniformOutput',false);
singles_count_g2_ti(k,3)=sum(cellfun('length',v2_g2_3));
singles_rate_g2_ti(k,3)=sum(cellfun('length',v2_g2_3))/(v_rate(3)*(60));

%% HOUR 4
%Find single fast and slow events
v2_g1_4=cellfun(@(equis1,equis2) single_hfo_get_sample(equis1,equis2),Mx_cortex_g1_ti_4,cohfos2_g1_4,'UniformOutput',false);
singles_count_g1_ti(k,4)=sum(cellfun('length',v2_g1_4));
singles_rate_g1_ti(k,4)=sum(cellfun('length',v2_g1_4))/(v_rate(4)*(60));


v2_g2_4=cellfun(@(equis1,equis2) single_hfo_get_sample(equis1,equis2),Mx_cortex_g2_ti_4,cohfos2_g2_4,'UniformOutput',false);
singles_count_g2_ti(k,4)=sum(cellfun('length',v2_g2_4));
singles_rate_g2_ti(k,4)=sum(cellfun('length',v2_g2_4))/(v_rate(4)*(60));
%xo

%%

%HPC COHFOS (Ripples coocurring with hfos)
cohf_mx_hpc=Mx_hpc(~cellfun('isempty',cohfos1));%Peak values cells where HPC cohfos were found.
cohf_sx_hpc=Sx_hpc(~cellfun('isempty',cohfos1));%Peak values cells where HPC cohfos were found.
cohf_ex_hpc=Ex_hpc(~cellfun('isempty',cohfos1));%Peak values cells where HPC cohfos were found.

Cohfos1=cohfos1(~cellfun('isempty',cohfos1));

%Locate sample per cohfos
coh_samp_hpc= cellfun(@(equis1,equis2) co_hfo_get_sample(equis1,equis2),cohf_mx_hpc,Cohfos1,'UniformOutput',false);

cohf_sx_hpc_val=cellfun(@(equis1,equis2) equis1(equis2),cohf_sx_hpc,coh_samp_hpc,'UniformOutput',false);
cohf_sx_hpc_val=[cohf_sx_hpc_val{:}];

cohf_mx_hpc_val=cellfun(@(equis1,equis2) equis1(equis2),cohf_mx_hpc,coh_samp_hpc,'UniformOutput',false);
cohf_mx_hpc_val=[cohf_mx_hpc_val{:}];

cohf_ex_hpc_val=cellfun(@(equis1,equis2) equis1(equis2),cohf_ex_hpc,coh_samp_hpc,'UniformOutput',false);
cohf_ex_hpc_val=[cohf_ex_hpc_val{:}];

%Duration of cohfos
cohf_hpc_dura=cohf_ex_hpc_val-cohf_sx_hpc_val;
cohf_hpc_dura=median(cohf_hpc_dura);
Cohf_hpc_dura(k)=cohf_hpc_dura;

%Cohfos ripples traces
Sig_hpc=sig_hpc(~cellfun('isempty',cohfos1));
Sig_hpc=cellfun(@(equis1,equis2) equis1(equis2),Sig_hpc,coh_samp_hpc,'UniformOutput',false);
Sig_hpc=[Sig_hpc{:}];


%Extract main features of coocuring ripples
[x,y,z,w,h,q,l,p]=hfo_specs(Sig_hpc,timeasleep,0,Rat,tr,fn);
fi_cohfo_hpc(k)=x;
fa_cohfo_hpc(k)=y;
amp_cohfo_hpc(k)=z;
count_cohfo_hpc(k)=w;
rate_cohfo_hpc(k)=h;
dura_cohfo_hpc(k)=q;
auc_cohfo_hpc(k)=l;
p2p_cohfo_hpc(k)=p;
%xo
%Single HPC ripples
v2=cellfun(@(equis1,equis2) single_hfo_get_sample(equis1,equis2),Mx_hpc,cohfos1,'UniformOutput',false);

Sig_hpc_single=cellfun(@(equis1,equis2) equis1(equis2),sig_hpc,v2,'UniformOutput',false);
Sig_hpc_single=[Sig_hpc_single{:}];

[single_mx_hpc_val,single_sx_hpc_val]=cellfun(@(equis1,equis2,equis3) single_hfos_mx(equis1,equis2,equis3),cohfos1,Mx_hpc,Sx_hpc,'UniformOutput',false);
single_mx_hpc_val=[single_mx_hpc_val{:}];
single_sx_hpc_val=[single_sx_hpc_val{:}];

[Sig_hpc_single_1]=find_single_events(Mx_ti_1,cohfos1_ti_1, sig_hpc);
[Sig_hpc_single_2]=find_single_events(Mx_ti_2,cohfos1_ti_2, sig_hpc);
[Sig_hpc_single_3]=find_single_events(Mx_ti_3,cohfos1_ti_3, sig_hpc);
[Sig_hpc_single_4]=find_single_events(Mx_ti_4,cohfos1_ti_4, sig_hpc);

%Compute main features of single ripples
[x,y,z,w,h,q,l,p]=hfo_specs(Sig_hpc_single,timeasleep,0,Rat,tr,fn);
fi_single_hpc(k)=x;
fa_single_hpc(k)=y;
amp_single_hpc(k)=z;
count_single_hpc(k)=w;
rate_single_hpc(k)=h;
dura_single_hpc(k)=q;
auc_single_hpc(k)=l;
p2p_single_hpc(k)=p;

count_single_hpc_perhour(k,1)=length(Sig_hpc_single_1);
count_single_hpc_perhour(k,2)=length(Sig_hpc_single_2);
count_single_hpc_perhour(k,3)=length(Sig_hpc_single_3);
count_single_hpc_perhour(k,4)=length(Sig_hpc_single_4);

rate_single_hpc_perhour(k,1)=length(Sig_hpc_single_1)/(v_rate(1)*60);
rate_single_hpc_perhour(k,2)=length(Sig_hpc_single_2)/(v_rate(2)*60);
rate_single_hpc_perhour(k,3)=length(Sig_hpc_single_3)/(v_rate(3)*60);
rate_single_hpc_perhour(k,4)=length(Sig_hpc_single_4)/(v_rate(4)*60);



%%%%
%Cortical COHFOS
cohf_mx_cortex=Mx_cortex(~cellfun('isempty',cohfos2));%Peak values cells where cortex cohfos were found.
cohf_sx_cortex=Sx_cortex(~cellfun('isempty',cohfos2));%Peak values cells where cortex cohfos were found.
cohf_ex_cortex=Ex_cortex(~cellfun('isempty',cohfos2));%Peak values cells where cortex cohfos were found.

Cohfos2=cohfos2(~cellfun('isempty',cohfos2));

%Locate sample per cohfos
coh_samp_cortex= cellfun(@(equis1,equis2) co_hfo_get_sample(equis1,equis2),cohf_mx_cortex,Cohfos2,'UniformOutput',false);
cohf_sx_cortex_val=cellfun(@(equis1,equis2) equis1(equis2),cohf_sx_cortex,coh_samp_cortex,'UniformOutput',false);
cohf_sx_cortex_val=[cohf_sx_cortex_val{:}];

cohf_mx_cortex_val=cellfun(@(equis1,equis2) equis1(equis2),cohf_mx_cortex,coh_samp_cortex,'UniformOutput',false);
cohf_mx_cortex_val=[cohf_mx_cortex_val{:}];

cohf_ex_cortex_val=cellfun(@(equis1,equis2) equis1(equis2),cohf_ex_cortex,coh_samp_cortex,'UniformOutput',false);
cohf_ex_cortex_val=[cohf_ex_cortex_val{:}];

%Duration of PPC cohfos
cohf_cortex_dura=cohf_ex_cortex_val-cohf_sx_cortex_val;
cohf_cortex_dura=median(cohf_cortex_dura);
Cohf_cortex_dura(k)=cohf_cortex_dura;

Sig_cortex=sig_cortex(~cellfun('isempty',cohfos2));
Sig_cortex=cellfun(@(equis1,equis2) equis1(equis2),Sig_cortex,coh_samp_cortex,'UniformOutput',false);
Sig_cortex=[Sig_cortex{:}];

%Compute features of PPC cohfos
[x,y,z,w,h,q,l,p]=hfo_specs(Sig_cortex,timeasleep,0,Rat,tr,fn);
fi_cohfo_cortex(k)=x;
fa_cohfo_cortex(k)=y;
amp_cohfo_cortex(k)=z;
count_cohfo_cortex(k)=w;
rate_cohfo_cortex(k)=h;
dura_cohfo_cortex(k)=q;
auc_cohfo_cortex(k)=l;
p2p_cohfo_cortex(k)=p;

%Single HFOs Cortex
v2=cellfun(@(equis1,equis2) single_hfo_get_sample(equis1,equis2),Mx_cortex,cohfos2,'UniformOutput',false);

Sig_cortex_single=cellfun(@(equis1,equis2) equis1(equis2),sig_cortex,v2,'UniformOutput',false);
Sig_cortex_single=[Sig_cortex_single{:}];

[single_mx_cortex_val,single_sx_cortex_val]=cellfun(@(equis1,equis2,equis3) single_hfos_mx(equis1,equis2,equis3),cohfos2,Mx_cortex,Sx_cortex,'UniformOutput',false);
single_mx_cortex_val=[single_mx_cortex_val{:}];
single_sx_cortex_val=[single_sx_cortex_val{:}];


[Sig_cortex_single_1]=find_single_events(Mx_cortex_ti_1,cohfos2_ti_1, sig_cortex);
[Sig_cortex_single_2]=find_single_events(Mx_cortex_ti_2,cohfos2_ti_2, sig_cortex);
[Sig_cortex_single_3]=find_single_events(Mx_cortex_ti_3,cohfos2_ti_3, sig_cortex);
[Sig_cortex_single_4]=find_single_events(Mx_cortex_ti_4,cohfos2_ti_4, sig_cortex);


%Compute features of single hfos
[x,y,z,w,h,q,l,p]=hfo_specs(Sig_cortex_single,timeasleep,0,Rat,tr,fn);
fi_single_cortex(k)=x;
fa_single_cortex(k)=y;
amp_single_cortex(k)=z;
count_single_cortex(k)=w;
rate_single_cortex(k)=h;
dura_single_cortex(k)=q;
auc_single_cortex(k)=l;
p2p_single_cortex(k)=p;

count_single_cortex_perhour(k,1)=length(Sig_cortex_single_1);
count_single_cortex_perhour(k,2)=length(Sig_cortex_single_2);
count_single_cortex_perhour(k,3)=length(Sig_cortex_single_3);
count_single_cortex_perhour(k,4)=length(Sig_cortex_single_4);

rate_single_cortex_perhour(k,1)=length(Sig_cortex_single_1)/(v_rate(1)*60);
rate_single_cortex_perhour(k,2)=length(Sig_cortex_single_2)/(v_rate(2)*60);
rate_single_cortex_perhour(k,3)=length(Sig_cortex_single_3)/(v_rate(3)*60);
rate_single_cortex_perhour(k,4)=length(Sig_cortex_single_4)/(v_rate(4)*60);

progress_bar(k,length(g),f)
    cd ..    
    end
xo        
%% Generate tables and save values into spreadsheets.

%%
%Cortex cohfos
    TT=table;
    TT.Variables=    [[{'Count'};{'Rate'};{'Duration'};{'Average Frequency'};{'Instantaneous Frequency'};{'Amplitude'}] num2cell([count_cohfo_cortex;rate_cohfo_cortex;dura_cohfo_cortex;fa_cohfo_cortex;fi_cohfo_cortex; amp_cohfo_cortex])];
    
    TT.Properties.VariableNames=['Metric';cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)].';
    
            writetable(TT,strcat(xx{1},'_',num2str(tr(2)),'_cohfos','.xls'),'Sheet',1,'Range','A2:L10')    

%Cortex singles
    TT=table;
    TT.Variables=    [[{'Count'};{'Rate'};{'Duration'};{'Average Frequency'};{'Instantaneous Frequency'};{'Amplitude'}] num2cell([count_single_cortex;rate_single_cortex;dura_single_cortex;fa_single_cortex;fi_single_cortex; amp_single_cortex])];
    
    TT.Properties.VariableNames=['Metric';cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)].';
    
            writetable(TT,strcat(xx{1},'_',num2str(tr(2)),'_singles','.xls'),'Sheet',1,'Range','A2:L10')    

%%
%HPC Multiplets
t1=repmat({'x'},[1 length(g)+2]);

for ll=1:length(multiplets)
    TT=table;
    eval(strcat('TT.Variables=    [','[','{' ,'''',multiplets{ll},'''','};','{' ,'''','x','''','};','{' ,'''','x','''','}',']'," ",'[','{' ,'''','Count','''','};','{' ,'''','Rate','''','};','{' ,'''','Duration','''','}',']',...
    ' num2cell([hfos_hpc_',multiplets{ll},';hfos_hpc_rate_',multiplets{ll},';hfos_hpc_duration_',multiplets{ll},'])];'))
    TT.Properties.VariableNames=['Event';'Metric';cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)].';

    tab.(multiplets{ll})=TT;
    if ll==1
        Tab=tab.(multiplets{ll});
    else
        Tab=[Tab;t1;tab.(multiplets{ll})];
    end

end

writetable(Tab,strcat('HPC','_',num2str(tr(1)),'_multiplets','.xls'),'Sheet',1,'Range','A1:Z50')


%%
%HPC ripples coocurrent with hfos (cohfos)
    TT=table;
    TT.Variables=    [[{'Count'};{'Rate'};{'Duration'};{'Average Frequency'};{'Instantaneous Frequency'};{'Amplitude'}] num2cell([count_cohfo_hpc;rate_cohfo_hpc;dura_cohfo_hpc;fa_cohfo_hpc;fi_cohfo_hpc; amp_cohfo_hpc])];
    
    TT.Properties.VariableNames=['Metric';cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)].';
    
            writetable(TT,strcat('HPC','_',num2str(tr(1)),'_',num2str(tr(2)),'_cohfos','.xls'),'Sheet',1,'Range','A2:L10')    

%HPC singles
    TT=table;
    TT.Variables=    [[{'Count'};{'Rate'};{'Duration'};{'Average Frequency'};{'Instantaneous Frequency'};{'Amplitude'}] num2cell([count_single_hpc;rate_single_hpc;dura_single_hpc;fa_single_hpc;fi_single_hpc; amp_single_hpc])];
    
    TT.Properties.VariableNames=['Metric';cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)].';
    
            writetable(TT,strcat('HPC','_',num2str(tr(1)),'_',num2str(tr(2)),'_singles','.xls'),'Sheet',1,'Range','A2:L10')    



%%  Cohfos
 
    TT=table;
    TT.Variables=    [[{'Count'};{'Rate'}] num2cell([cohfos_count;cohfos_rate;])];
    TT.Properties.VariableNames=['Metric';g];    
    writetable(TT,strcat('coHFOs_',num2str(tr(2)),'.xls'),'Sheet',1,'Range','A2:L6')    
    
%% Slow and fast hfos 

%coocurrent HFOs
    TT=table;
    TT.Variables=    [[{'Slower'};{'x'};{'Faster'};{'x'}] [{'Count'};{'Rate'};{'Count'};{'Rate'}] num2cell([cohfos_count_g1;cohfos_rate_g1;cohfos_count_g2;cohfos_rate_g2])];

    TT.Properties.VariableNames=['Events';'Metric';g];


writetable(TT,strcat('slower_faster_cohfos_',num2str(tr(1)),'_',num2str(tr(2)),'.xls'),'Sheet',1,'Range','A1:Z50')

%Single HFOs
    TT=table;
    TT.Variables=    [[{'Slower'};{'x'};{'Faster'};{'x'}] [{'Count'};{'Rate'};{'Count'};{'Rate'}] num2cell([singles_count_g1;singles_rate_g1;singles_count_g2;singles_rate_g2])];

    TT.Properties.VariableNames=['Events';'Metric';g];


writetable(TT,strcat('slower_faster_singles_',num2str(tr(1)),'_',num2str(tr(2)),'.xls'),'Sheet',1,'Range','A1:Z50')

%% New Section: Counts and rates per hour.            
%%  Cohfos per hour
 
%Counts
    TT=table;
    TT.Variables=    [[{'Hour1'};{'Hour2'};{'Hour3'};{'Hour4'};{'Total'}] num2cell([cohfos_count_perhour.';sum(cohfos_count_perhour.')])];
    TT.Properties.VariableNames=[{'Hour'};g];
    writetable(TT,strcat('cohfos_counts_perhour_',num2str(tr(1)),'_',num2str(tr(2)),'.xls'),'Sheet',1,'Range','A1:Z50')
    
%Rates
    TT=table;
    TT.Variables=    [[{'Hour1'};{'Hour2'};{'Hour3'};{'Hour4'}] num2cell([cohfos_rate_perhour.'])];
    TT.Properties.VariableNames=[{'Hour'};g];
    writetable(TT,strcat('cohfos_rate_perhour_',num2str(tr(1)),'_',num2str(tr(2)),'.xls'),'Sheet',1,'Range','A1:Z50')
%% Singles per hour
%HPC 

%Counts
    TT=table;
    TT.Variables=    [[{'Hour1'};{'Hour2'};{'Hour3'};{'Hour4'};{'Total'}] num2cell([count_single_hpc_perhour.';sum(count_single_hpc_perhour.')])];
    TT.Properties.VariableNames=[{'Hour'};g];
    writetable(TT,strcat('singles_hpc_counts_perhour_',num2str(tr(1)),'_',num2str(tr(2)),'.xls'),'Sheet',1,'Range','A1:Z50')
    
%Rates
    TT=table;
    TT.Variables=    [[{'Hour1'};{'Hour2'};{'Hour3'};{'Hour4'}] num2cell([rate_single_hpc_perhour.'])];
    TT.Properties.VariableNames=[{'Hour'};g];
    writetable(TT,strcat('singles_hpc_rate_perhour_',num2str(tr(1)),'_',num2str(tr(2)),'.xls'),'Sheet',1,'Range','A1:Z50')

%CORTEX
%Counts
    TT=table;
    TT.Variables=    [[{'Hour1'};{'Hour2'};{'Hour3'};{'Hour4'};{'Total'}] num2cell([count_single_cortex_perhour.';sum(count_single_cortex_perhour.')])];
    TT.Properties.VariableNames=[{'Hour'};g];
    writetable(TT,strcat('singles_cortex_counts_perhour_',num2str(tr(1)),'_',num2str(tr(2)),'.xls'),'Sheet',1,'Range','A1:Z50')
    
%Rates
    TT=table;
    TT.Variables=    [[{'Hour1'};{'Hour2'};{'Hour3'};{'Hour4'}] num2cell([rate_single_cortex_perhour.'])];
    TT.Properties.VariableNames=[{'Hour'};g];
    writetable(TT,strcat('singles_cortex_rate_perhour_',num2str(tr(1)),'_',num2str(tr(2)),'.xls'),'Sheet',1,'Range','A1:Z50')


%% Slow and fast hfos 

%coocurrent HFOs

%Slow Counts
    TT=table;
    TT.Variables=    [[{'Hour1'};{'Hour2'};{'Hour3'};{'Hour4'};{'Total'}] num2cell([cohfos_count_g1_ti.';sum(cohfos_count_g1_ti.')])];
    TT.Properties.VariableNames=[{'Hour'};g];
    writetable(TT,strcat('slow_cohfos_counts_perhour_',num2str(tr(1)),'_',num2str(tr(2)),'.xls'),'Sheet',1,'Range','A1:Z50')

%Fast Counts
    TT=table;
    TT.Variables=    [[{'Hour1'};{'Hour2'};{'Hour3'};{'Hour4'};{'Total'}] num2cell([cohfos_count_g2_ti.';sum(cohfos_count_g2_ti.')])];
    TT.Properties.VariableNames=[{'Hour'};g];
    writetable(TT,strcat('fast_cohfos_counts_perhour_',num2str(tr(1)),'_',num2str(tr(2)),'.xls'),'Sheet',1,'Range','A1:Z50')

%Slow Rates
    TT=table;
    TT.Variables=    [[{'Hour1'};{'Hour2'};{'Hour3'};{'Hour4'}] num2cell([cohfos_rate_g1_ti.'])];
    TT.Properties.VariableNames=[{'Hour'};g];
    writetable(TT,strcat('slow_cohfos_rates_perhour_',num2str(tr(1)),'_',num2str(tr(2)),'.xls'),'Sheet',1,'Range','A1:Z50')

%Fast Rates
    TT=table;
    TT.Variables=    [[{'Hour1'};{'Hour2'};{'Hour3'};{'Hour4'}] num2cell([cohfos_rate_g2_ti.'])];
    TT.Properties.VariableNames=[{'Hour'};g];
    writetable(TT,strcat('fast_cohfos_rates_perhour_',num2str(tr(1)),'_',num2str(tr(2)),'.xls'),'Sheet',1,'Range','A1:Z50')
    
%% 


%Single HFOs Counts G1
TT=table;
TT.Variables=    [[{'Hour1'};{'Hour2'};{'Hour3'};{'Hour4'};{'Total'}] num2cell([singles_count_g1_ti.';sum(singles_count_g1_ti.')])];
TT.Properties.VariableNames=[{'Hour'};g];
writetable(TT,strcat('slow_single_counts_perhour_',num2str(tr(1)),'_',num2str(tr(2)),'.xls'),'Sheet',1,'Range','A1:Z50')


%Single HFOs Counts G2    
TT=table;
TT.Variables=    [[{'Hour1'};{'Hour2'};{'Hour3'};{'Hour4'};{'Total'}] num2cell([singles_count_g2_ti.';sum(singles_count_g2_ti.')])];
TT.Properties.VariableNames=[{'Hour'};g];
writetable(TT,strcat('fast_single_counts_perhour_',num2str(tr(1)),'_',num2str(tr(2)),'.xls'),'Sheet',1,'Range','A1:Z50')


%Single HFOs Rates G1    
TT=table;
TT.Variables=    [[{'Hour1'};{'Hour2'};{'Hour3'};{'Hour4'}] num2cell([singles_rate_g1_ti.'])];
TT.Properties.VariableNames=[{'Hour'};g];
writetable(TT,strcat('slow_single_rates_perhour_',num2str(tr(1)),'_',num2str(tr(2)),'.xls'),'Sheet',1,'Range','A1:Z50')


%Single HFOs Rates G2    
TT=table;
TT.Variables=    [[{'Hour1'};{'Hour2'};{'Hour3'};{'Hour4'}] num2cell([singles_rate_g2_ti.'])];
TT.Properties.VariableNames=[{'Hour'};g];
writetable(TT,strcat('fast_single_rates_perhour_',num2str(tr(1)),'_',num2str(tr(2)),'.xls'),'Sheet',1,'Range','A1:Z50')

%% 30_minute counts

for i=1:length(g)
    tet=length(CR_30_count{1,i});
    if tet<10
        %     C_temp{i} = CR_30_count{1,i};
        CR_30_count{1,i}(1,tet+1:10)= zeros(1,10-tet);
    elseif tet>10
        xo % if this happens, modify the code to include the excess time
    end
end

for i=1:length(g)
    tet=length(HR_30_count{1,i});
    if tet<10
        %     H_temp{i} = HR_30_count{1,i};
        HR_30_count{1,i}(1,tet+1:10)= zeros(1,10-tet);
    elseif tet>10
        xo % if this happens, modify the code to include the excess time
    end
end
    temp = [num2cell(CR_30_count{1,1}');num2cell(HR_30_count{1,1}')];
for k=2:length(g)
        tt = [num2cell(CR_30_count{1,k}');num2cell(HR_30_count{1,k}')];
        temp = [temp,tt];
end
    TT=table;
    TT.Variables=    [[{'Cortical Ripples count'};{''};{''};{''};...
        {''};{''};{''};{''};{''};{''};...
        {'Hippocampal Ripples count'};{''};{''};{''};...
        {''};{''};{''};{''};{''};{''}]...
        [{'30 minutes '};{'1 hour'};{'1.5 hours'};{'2 hours'};...
        {'2.5 hours'};{'3 hours'};{'3.5 hours '};{'4 hours'};...
        {'4.5 hours'};{'5 hours'};...
        {'30 minutes '};{'1 hour'};{'1.5 hours'};{'2 hours'};...
        {'2.5 hours'};{'3 hours'};{'3.5 hours '};{'4 hours'};...
        {'4.5 hours'};{'5 hours'}],...
        [temp]];
    
    TT.Properties.VariableNames=['Ripple type';'Duration';g];

    writetable(TT,strcat('30_min_count_',num2str(cell2mat(answer)),'.xls'),'Sheet',1,'Range','A2:L24')

%% NREM duration
    for i=1:length(g)
        tet=length(v_rate_CR{1,i});
        if tet<10
                v_rate_CR_t{i} = v_rate_CR{1,i};
            v_rate_CR{1,i}(tet+1:10,1)= zeros(10-tet,1);
        elseif tet>10
            xo % if this happens, modify the code to include the excess time
        end
    end

    temp = [num2cell(v_rate_CR{1,1})]; %;num2cell(v_rate_HR{1,1});num2cell(v_rate_HR{1,k})
    for k=2:length(g)
        tt = [num2cell(v_rate_CR{1,k})];
        temp = [temp,tt];
    end
    TT=table;
    TT.Variables=    [[{'NREM duration'};{''};{''};{''};...
        {''};{''};{''};{''};{''};{''}]...
        [{'30 minutes '};{'1 hour'};{'1.5 hours'};{'2 hours'};...
        {'2.5 hours'};{'3 hours'};{'3.5 hours '};{'4 hours'};...
        {'4.5 hours'};{'5 hours'}]...
        [temp]];
    TT.Properties.VariableNames=['NREM';'Duration';g];

    writetable(TT,strcat('30_min_NREM_',num2str(cell2mat(answer)),'.xls'),'Sheet',1,'Range','A2:L20')

return
