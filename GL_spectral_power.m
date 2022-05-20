%GL_spectral_power

% Main code for computing event-centered spectral power.
% Requires 'load_me_first.mat' loaded first. 

%% Find location
close all
dname=uigetdir([],'Select folder with Matlab data containing all rats.');
cd(dname)

% cd('/home/adrian/Documents/Plusmaze_downsampled')

%%
%Select rat number
opts.Resize = 'on';
opts.WindowStyle = 'modal';
opts.Interpreter = 'tex';
prompt=strcat('\bf Select a rat#. Options:','{ }',num2str(rats));
answer = inputdlg(prompt,'Input',[2 30],{''},opts);
Rat=str2num(answer{1});
cd(num2str(Rat))
tr=getfield(T,strcat('Rat',num2str(Rat)));%Thresholds 
xx={'PAR'};% Posterior Parietal cortex.
%% Get folder names
labelconditions=getfolder;
labelconditions=labelconditions.';
g=labelconditions;
multiplets=[{'singlets'} {'doublets'} {'triplets'} {'quatruplets'} {'pentuplets'} {'sextuplets'} {'septuplets'} {'octuplets'} {'nonuplets'}];
iii=1;
fn=1000; %Sampling frequency after downsampling.
%% Display and select conditions and sessions to analyze
    %Center figure.
    f=figure();
    movegui(gcf,'center');

    %Checkboxes
    Boxcheck = cell(1,4);
    for h1=1:length(labelconditions)
    boxcheck = uicontrol(f,'Style','checkbox','String',labelconditions{h1},'Position',[10 f.Position(4)-30*h1 400 20]);
    boxcheck.FontSize=11;
    boxcheck.Value=1;
    Boxcheck{h1}=boxcheck;   
    end

    set(f, 'NumberTitle', 'off', ...
        'Name', 'Select conditions and sessions');

    %Push button
    c = uicontrol;
    c.String = 'Continue';
    c.FontSize=10;
    c.Position=[f.Position(1)/3.5 c.Position(2)-10 f.Position(3)/2 c.Position(4)];

    %Callback
    c.Callback='uiresume(gcbf)';
    uiwait(gcf); 
    boxch=cellfun(@(x) get(x,'Value'),Boxcheck);
    clear Boxcheck
    close(f);
g={g{logical(boxch)}};    

if sum(cell2mat(cellfun(@(equis1) contains(equis1,'nl'),g,'UniformOutput',false)))==1
g=g([find(cell2mat(cellfun(@(equis1) contains(equis1,labelconditions2{1}),g,'UniformOutput',false)))...
 find(cell2mat(cellfun(@(equis1) contains(equis1,labelconditions2{2}),g,'UniformOutput',false)))...
 find(cell2mat(cellfun(@(equis1) contains(equis1,labelconditions2{3}),g,'UniformOutput',false)))...
 find(cell2mat(cellfun(@(equis1) contains(equis1,labelconditions2{4}),g,'UniformOutput',false)))]);

else
    error('Name issue')
end

%%
f=waitbar(0,'Please wait...');
    for k=1:length(g)
        cd(g{k})
% Load signals
CORTEX=dir(strcat('*',xx{1},'*.mat'));
if isempty(CORTEX)
    g=g(~contains(g,g{k}));s
    cd ..
    progress_bar(k,length(g),f)
    break
end

%PARIETAL
CORTEX=CORTEX.name;
CORTEX=load(CORTEX);
CORTEX=getfield(CORTEX,xx{1});
CORTEX=CORTEX.*(0.195); % Open Ephys BitVolt factor

%HPC
HPC=dir(strcat('*','HPC','*.mat'));
HPC=HPC.name;
HPC=load(HPC);
HPC=getfield(HPC,'HPC');
HPC=HPC.*(0.195);

%PFC
PFC=dir(strcat('*','PFC','*.mat'));
PFC=PFC.name;
PFC=load(PFC);
PFC=getfield(PFC,'PFC');
PFC=PFC.*(0.195);


%Load sleep stages
A = dir('*states*.mat');
A={A.name};

if  ~isempty(A)
       cellfun(@load,A);
else
      error('No Scoring found')    
end
 
 ro=150; %Select length of window in miliseconds
 %Find HFOss
[Mx_cortex,timeasleep,sig_cortex,Ex_cortex,Sx_cortex,...
  p_cortex,q_cortex,cont_cortex,sig_pq_cortex ...
  ]=gui_findripples_spec(CORTEX,states,xx,tr,PFC,HPC,fn,ro);
%
%p: raw signal
%q: bandpassed signal in ripple band.

%Get traces of events
si=sig_cortex(~cellfun('isempty',sig_cortex));
si=[si{:}];

si_pq=sig_pq_cortex(~cellfun('isempty',sig_pq_cortex));
si_pq=[si_pq{:}];

%Event-centered traces in the Fieldtrip format
Q_cortex=q_cortex(~cellfun('isempty',q_cortex));
Q_cortex=[Q_cortex{:}];

%Find slow and fast hfos.
print_hist=0;
[~,~,~,~,~,~,~,~,si_mixed,~]=hfo_specs(si,timeasleep,print_hist,Rat,tr);
%Find empty cells
void_index=find(cellfun('isempty',Q_cortex));

%All par HFOS splitted in slow and fast.
Q_cortex_g1=Q_cortex(si_mixed.i1(~ismember(si_mixed.i1,void_index)));
Q_cortex_g2=Q_cortex(si_mixed.i2(~ismember(si_mixed.i2,void_index)));

%% HPC     

%Find ripples and ripple-centered traces
[Mx_hpc,timeasleep,sig_hpc,Ex_hpc,Sx_hpc,...
  p_hpc,q_hpc,cont_hpc ...
]=gui_findripples_spec(HPC,states,{'HPC'},tr,PFC,CORTEX,fn,ro);

si=sig_hpc(~cellfun('isempty',sig_hpc));
si=[si{:}];

%% Coocurent hfos
[cohfos1,cohfos2]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),Mx_hpc,Mx_cortex,'UniformOutput',false);
%cohfos1: HPC timestamps.
%cohfos2: Cortex timestamps.

%%

%HPC COHFOS
cohf_mx_hpc=Mx_hpc(~cellfun('isempty',cohfos1));%Peak values cells where HPC cohfos were found.
cohf_sx_hpc=Sx_hpc(~cellfun('isempty',cohfos1));%Start values cells where HPC cohfos were found.
cohf_ex_hpc=Ex_hpc(~cellfun('isempty',cohfos1));%End values cells where HPC cohfos were found.

Cohfos1=cohfos1(~cellfun('isempty',cohfos1));

%Locate sample per cohfos
coh_samp_hpc= cellfun(@(equis1,equis2) co_hfo_get_sample(equis1,equis2),cohf_mx_hpc,Cohfos1,'UniformOutput',false);

cohf_sx_hpc_val=cellfun(@(equis1,equis2) equis1(equis2),cohf_sx_hpc,coh_samp_hpc,'UniformOutput',false);
cohf_sx_hpc_val=[cohf_sx_hpc_val{:}];

cohf_mx_hpc_val=cellfun(@(equis1,equis2) equis1(equis2),cohf_mx_hpc,coh_samp_hpc,'UniformOutput',false);
cohf_mx_hpc_val=[cohf_mx_hpc_val{:}];

cohf_ex_hpc_val=cellfun(@(equis1,equis2) equis1(equis2),cohf_ex_hpc,coh_samp_hpc,'UniformOutput',false);
cohf_ex_hpc_val=[cohf_ex_hpc_val{:}];

cohf_hpc_dura=cohf_ex_hpc_val-cohf_sx_hpc_val;
cohf_hpc_dura=median(cohf_hpc_dura);
Cohf_hpc_dura(k)=cohf_hpc_dura;

Sig_hpc=sig_hpc(~cellfun('isempty',cohfos1));
Sig_hpc=cellfun(@(equis1,equis2) equis1(equis2),Sig_hpc,coh_samp_hpc,'UniformOutput',false);
Sig_hpc=[Sig_hpc{:}];

%COHFOS windows
p_cohfos_hpc=p_hpc(~cellfun('isempty',cohfos1));
p_cohfos_hpc=cellfun(@(equis1,equis2) equis1(equis2),p_cohfos_hpc,coh_samp_hpc,'UniformOutput',false);
p_cohfos_hpc=[p_cohfos_hpc{:}];
q_cohfos_hpc=q_hpc(~cellfun('isempty',cohfos1));
q_cohfos_hpc=cellfun(@(equis1,equis2) equis1(equis2),q_cohfos_hpc,coh_samp_hpc,'UniformOutput',false);
q_cohfos_hpc=[q_cohfos_hpc{:}];
p_cohfos_hpc=p_cohfos_hpc(~cellfun('isempty',p_cohfos_hpc));
q_cohfos_hpc=q_cohfos_hpc(~cellfun('isempty',q_cohfos_hpc));


%Single events
v2=cellfun(@(equis1,equis2) single_hfo_get_sample(equis1,equis2),Mx_hpc,cohfos1,'UniformOutput',false);

Sig_hpc_single=cellfun(@(equis1,equis2) equis1(equis2),sig_hpc,v2,'UniformOutput',false);
Sig_hpc_single=[Sig_hpc_single{:}];

%single windows
p_single_hpc=cellfun(@(equis1,equis2) equis1(equis2),p_hpc,v2,'UniformOutput',false);
p_single_hpc=[p_single_hpc{:}];
q_single_hpc=cellfun(@(equis1,equis2) equis1(equis2),q_hpc,v2,'UniformOutput',false);
q_single_hpc=[q_single_hpc{:}];
p_single_hpc=p_single_hpc(~cellfun('isempty',p_single_hpc));
q_single_hpc=q_single_hpc(~cellfun('isempty',q_single_hpc));


[single_mx_hpc_val,single_sx_hpc_val]=cellfun(@(equis1,equis2,equis3) single_hfos_mx(equis1,equis2,equis3),cohfos1,Mx_hpc,Sx_hpc,'UniformOutput',false);
single_mx_hpc_val=[single_mx_hpc_val{:}];
single_sx_hpc_val=[single_sx_hpc_val{:}];


%%%%
%Cortical COHFOS
cohf_mx_cortex=Mx_cortex(~cellfun('isempty',cohfos2));%Peak values cells where cortex cohfos were found.
cohf_sx_cortex=Sx_cortex(~cellfun('isempty',cohfos2));%Start values cells where cortex cohfos were found.
cohf_ex_cortex=Ex_cortex(~cellfun('isempty',cohfos2));%End values cells where cortex cohfos were found.

Cohfos2=cohfos2(~cellfun('isempty',cohfos2));

%Locate sample per cohfos
coh_samp_cortex= cellfun(@(equis1,equis2) co_hfo_get_sample(equis1,equis2),cohf_mx_cortex,Cohfos2,'UniformOutput',false);
cohf_sx_cortex_val=cellfun(@(equis1,equis2) equis1(equis2),cohf_sx_cortex,coh_samp_cortex,'UniformOutput',false);
cohf_sx_cortex_val=[cohf_sx_cortex_val{:}];

cohf_mx_cortex_val=cellfun(@(equis1,equis2) equis1(equis2),cohf_mx_cortex,coh_samp_cortex,'UniformOutput',false);
cohf_mx_cortex_val=[cohf_mx_cortex_val{:}];

cohf_ex_cortex_val=cellfun(@(equis1,equis2) equis1(equis2),cohf_ex_cortex,coh_samp_cortex,'UniformOutput',false);
cohf_ex_cortex_val=[cohf_ex_cortex_val{:}];
cohf_cortex_dura=cohf_ex_cortex_val-cohf_sx_cortex_val;
cohf_cortex_dura=median(cohf_cortex_dura);
Cohf_cortex_dura(k)=cohf_cortex_dura;

Sig_cortex=sig_cortex(~cellfun('isempty',cohfos2));
Sig_cortex=cellfun(@(equis1,equis2) equis1(equis2),Sig_cortex,coh_samp_cortex,'UniformOutput',false);
Sig_cortex=[Sig_cortex{:}];

%COHFOS windows
p_cohfos_cortex=p_cortex(~cellfun('isempty',cohfos2));
p_cohfos_cortex=cellfun(@(equis1,equis2) equis1(equis2),p_cohfos_cortex,coh_samp_cortex,'UniformOutput',false);
p_cohfos_cortex=[p_cohfos_cortex{:}];
q_cohfos_cortex=q_cortex(~cellfun('isempty',cohfos2));
q_cohfos_cortex=cellfun(@(equis1,equis2) equis1(equis2),q_cohfos_cortex,coh_samp_cortex,'UniformOutput',false);
q_cohfos_cortex=[q_cohfos_cortex{:}];
p_cohfos_cortex=p_cohfos_cortex(~cellfun('isempty',p_cohfos_cortex));
q_cohfos_cortex=q_cohfos_cortex(~cellfun('isempty',q_cohfos_cortex));

%Single events
v2=cellfun(@(equis1,equis2) single_hfo_get_sample(equis1,equis2),Mx_cortex,cohfos2,'UniformOutput',false);

Sig_cortex_single=cellfun(@(equis1,equis2) equis1(equis2),sig_cortex,v2,'UniformOutput',false);
Sig_cortex_single=[Sig_cortex_single{:}];

%Single cortex windows
p_single_cortex=cellfun(@(equis1,equis2) equis1(equis2),p_cortex,v2,'UniformOutput',false);
p_single_cortex=[p_single_cortex{:}];
q_single_cortex=cellfun(@(equis1,equis2) equis1(equis2),q_cortex,v2,'UniformOutput',false);
q_single_cortex=[q_single_cortex{:}];
p_single_cortex=p_single_cortex(~cellfun('isempty',p_single_cortex));
q_single_cortex=q_single_cortex(~cellfun('isempty',q_single_cortex));


[single_mx_cortex_val,single_sx_cortex_val]=cellfun(@(equis1,equis2,equis3) single_hfos_mx(equis1,equis2,equis3),cohfos2,Mx_cortex,Sx_cortex,'UniformOutput',false);
single_mx_cortex_val=[single_mx_cortex_val{:}];
single_sx_cortex_val=[single_sx_cortex_val{:}];


 p_hpc=p_hpc(~cellfun('isempty',p_hpc));
 p_hpc=[p_hpc{:}];
 q_hpc=q_hpc(~cellfun('isempty',q_hpc));
 q_hpc=[q_hpc{:}];
 p_hpc=p_hpc(~cellfun('isempty',p_hpc));
 q_hpc=q_hpc(~cellfun('isempty',q_hpc));
 

%% Mixed distribution (Average freq) coHFOs
%Initialize variables
Mx_cortex_g1=Mx_cortex;
Mx_cortex_g2=Mx_cortex;
Ex_cortex_g1=Ex_cortex;
Ex_cortex_g2=Ex_cortex;
Sx_cortex_g1=Sx_cortex;
Sx_cortex_g2=Sx_cortex;

row=si_mixed.i1; %Indices for slow hfos
cont=0; %Counter
for ll=1:length(Mx_cortex)

    if ~isempty(Mx_cortex{ll})

        for lll=1:length(Mx_cortex{ll})
            cont=cont+1;

            if ~ismember(cont,row)
                Mx_cortex_g1{ll}(lll)=NaN;
                Ex_cortex_g1{ll}(lll)=NaN;
                Sx_cortex_g1{ll}(lll)=NaN;
                
            else
                Mx_cortex_g2{ll}(lll)=NaN;
                Ex_cortex_g2{ll}(lll)=NaN;
                Sx_cortex_g2{ll}(lll)=NaN;

            end

        end
         Mx_cortex_g1{ll}=Mx_cortex_g1{ll}(~isnan(Mx_cortex_g1{ll}));
         Mx_cortex_g2{ll}=Mx_cortex_g2{ll}(~isnan(Mx_cortex_g2{ll}));

         Ex_cortex_g1{ll}=Ex_cortex_g1{ll}(~isnan(Ex_cortex_g1{ll}));
         Ex_cortex_g2{ll}=Ex_cortex_g2{ll}(~isnan(Ex_cortex_g2{ll}));
         Sx_cortex_g1{ll}=Sx_cortex_g1{ll}(~isnan(Sx_cortex_g1{ll}));
         Sx_cortex_g2{ll}=Sx_cortex_g2{ll}(~isnan(Sx_cortex_g2{ll}));
         
         
    end

end
%Slow and fast cohfos
[~,cohfos2_g1]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),Mx_hpc,Mx_cortex_g1,'UniformOutput',false);
[~,cohfos2_g2]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),Mx_hpc,Mx_cortex_g2,'UniformOutput',false);

%Single HFOs
v2_g1=cellfun(@(equis1,equis2) single_hfo_get_sample(equis1,equis2),Mx_cortex_g1,cohfos2_g1,'UniformOutput',false);
v2_g2=cellfun(@(equis1,equis2) single_hfo_get_sample(equis1,equis2),Mx_cortex_g2,cohfos2_g2,'UniformOutput',false);


[p_cohfos_cortex_g1,q_cohfos_cortex_g1,p_single_cortex_g1,q_single_cortex_g1]=get_window_slowfast(Mx_cortex,Sx_cortex,Ex_cortex, cohfos2_g1,p_cortex,q_cortex,v2_g1);
[p_cohfos_cortex_g2,q_cohfos_cortex_g2,p_single_cortex_g2,q_single_cortex_g2]=get_window_slowfast(Mx_cortex,Sx_cortex,Ex_cortex, cohfos2_g2,p_cortex,q_cortex,v2_g2);

%%
  p_cortex=p_cortex(~cellfun('isempty',p_cortex));
 p_cortex=[p_cortex{:}];
 q_cortex=q_cortex(~cellfun('isempty',q_cortex));
 q_cortex=[q_cortex{:}];
 p_cortex=p_cortex(~cellfun('isempty',p_cortex));
 q_cortex=q_cortex(~cellfun('isempty',q_cortex));

%%
     P.(strrep(labelconditions2{k},'-','_')).(label1{1})={p_cohfos_hpc;p_single_hpc;p_hpc};
     P.(strrep(labelconditions2{k},'-','_')).(label1{3})={p_cohfos_cortex;p_single_cortex;p_cortex};
    
     Q.(strrep(labelconditions2{k},'-','_')).(label1{1})={q_cohfos_hpc;q_single_hpc;q_hpc};
     Q.(strrep(labelconditions2{k},'-','_')).(label1{3})={q_cohfos_cortex;q_single_cortex;q_cortex};
     
     %Slow/Fast-centered traces in the Fieldtrip format.
     %P: raw signal
     %Q: bandpassed signal in ripple band
     
     SP.(strrep(labelconditions2{k},'-','_')).(label1{3})={p_cohfos_cortex_g1;p_single_cortex_g1};
     SQ.(strrep(labelconditions2{k},'-','_')).(label1{3})={q_cohfos_cortex_g1;q_single_cortex_g1};
     
     FP.(strrep(labelconditions2{k},'-','_')).(label1{3})={p_cohfos_cortex_g2;p_single_cortex_g2};
     FQ.(strrep(labelconditions2{k},'-','_')).(label1{3})={q_cohfos_cortex_g2;q_single_cortex_g2};

progress_bar(k,length(g),f)
    cd ..    
    end
%Stop code here and save values in excel sheets.    
%xo

%% Find values for slow/fast PAR hfos
win_size=50; %Window size in miliseconds
%Different binary options:
random_hfo=0; %If selecting a random hfo.
rand_first_run=0; %If you run the random selection for the first time.
same_nr_types=0; %Same N number of events across types.

if same_nr_types==0
    N=[];
end

%SLOW
%PAR COHFOSs
s=1; % s=1 for cooccurrent HFOs, s=2 for single events.
w=3; %Index from 'label1'. w=1 for hpc centered, w=3 for par centered.
%Find mean power value of the window for different frequency bands.
[values_spec]=getval_spectra(SP,SQ,labelconditions2,label1,s,w,win_size,same_nr_types,N,random_hfo,rand_first_run,tr);
TT1=table;
TT1.Variables=    [[{'100-250Hz'};{'100-150Hz'};{'150-200Hz'};{'200-250Hz'}] num2cell([values_spec.nl(:,1) values_spec.plusmaze(:,1) values_spec.novelty(:,1) values_spec.for(:,1) values_spec.nl(:,2) values_spec.plusmaze(:,2) values_spec.novelty(:,2) values_spec.for(:,2) values_spec.nl(:,3) values_spec.plusmaze(:,3) values_spec.novelty(:,3) values_spec.for(:,3)])];
TT1.Properties.VariableNames=[{'Range'};{'HPC Baseline'};{'HPC Plusmaze'};{'HPC Novelty'};{'HPC Foraging'};{'PFC Baseline'};{'PFC Plusmaze'};{'PFC Novelty'};{'PFC Foraging'};{'PAR Baseline'};{'PAR Plusmaze'};{'PAR Novelty'};{'PAR Foraging'}];    

%SLOW 
%PAR singles
s=2;
w=3;
[values_spec]=getval_spectra(SP,SQ,labelconditions2,label1,s,w,win_size,same_nr_types,N,random_hfo,rand_first_run,tr);
TT2=table;
TT2.Variables=    [[{'100-250Hz'};{'100-150Hz'};{'150-200Hz'};{'200-250Hz'}] num2cell([values_spec.nl(:,1) values_spec.plusmaze(:,1) values_spec.novelty(:,1) values_spec.for(:,1) values_spec.nl(:,2) values_spec.plusmaze(:,2) values_spec.novelty(:,2) values_spec.for(:,2) values_spec.nl(:,3) values_spec.plusmaze(:,3) values_spec.novelty(:,3) values_spec.for(:,3)])];
TT2.Properties.VariableNames=[{'Range'};{'HPC Baseline'};{'HPC Plusmaze'};{'HPC Novelty'};{'HPC Foraging'};{'PFC Baseline'};{'PFC Plusmaze'};{'PFC Novelty'};{'PFC Foraging'};{'PAR Baseline'};{'PAR Plusmaze'};{'PAR Novelty'};{'PAR Foraging'}];    
 
%FAST
%PAR COHFOS
s=1;
w=3;
[values_spec,n3]=getval_spectra(FP,FQ,labelconditions2,label1,s,w,win_size,same_nr_types,N,random_hfo,rand_first_run,tr);
TT3=table;
TT3.Variables=    [[{'100-250Hz'};{'100-150Hz'};{'150-200Hz'};{'200-250Hz'}] num2cell([values_spec.nl(:,1) values_spec.plusmaze(:,1) values_spec.novelty(:,1) values_spec.for(:,1) values_spec.nl(:,2) values_spec.plusmaze(:,2) values_spec.novelty(:,2) values_spec.for(:,2) values_spec.nl(:,3) values_spec.plusmaze(:,3) values_spec.novelty(:,3) values_spec.for(:,3)])];
TT3.Properties.VariableNames=[{'Range'};{'HPC Baseline'};{'HPC Plusmaze'};{'HPC Novelty'};{'HPC Foraging'};{'PFC Baseline'};{'PFC Plusmaze'};{'PFC Novelty'};{'PFC Foraging'};{'PAR Baseline'};{'PAR Plusmaze'};{'PAR Novelty'};{'PAR Foraging'}];    

%FAST 
%PAR singles
s=2;
w=3;
[values_spec,n4]=getval_spectra(FP,FQ,labelconditions2,label1,s,w,win_size,same_nr_types,N,random_hfo,rand_first_run,tr);
TT4=table;
TT4.Variables=    [[{'100-250Hz'};{'100-150Hz'};{'150-200Hz'};{'200-250Hz'}] num2cell([values_spec.nl(:,1) values_spec.plusmaze(:,1) values_spec.novelty(:,1) values_spec.for(:,1) values_spec.nl(:,2) values_spec.plusmaze(:,2) values_spec.novelty(:,2) values_spec.for(:,2) values_spec.nl(:,3) values_spec.plusmaze(:,3) values_spec.novelty(:,3) values_spec.for(:,3)])];
TT4.Properties.VariableNames=[{'Range'};{'HPC Baseline'};{'HPC Plusmaze'};{'HPC Novelty'};{'HPC Foraging'};{'PFC Baseline'};{'PFC Plusmaze'};{'PFC Novelty'};{'PFC Foraging'};{'PAR Baseline'};{'PAR Plusmaze'};{'PAR Novelty'};{'PAR Foraging'}];    

t1=repmat({'x'},[1 13]);

tab=[TT1;t1;TT2;t1;TT3;t1;TT4];

if win_size== 25
writetable(tab,strcat('spec_SLOWFAST_values_25_rat_', num2str(Rat),'_' ,num2str(tr(2)),'.xls'),'Sheet',1,'Range','A2:Z50')  
end

if win_size== 50
writetable(tab,strcat('spec_SLOWFAST_values_rat_', num2str(Rat),'_' ,num2str(tr(2)),'.xls'),'Sheet',1,'Range','A2:Z50')  
end

%% 
% Using ALL events.
win_size=50;
%HPC COHFOS
s=1;
w=1;
[values_spec]=getval_spectra_All(P,Q,labelconditions2,label1,s,w,win_size);
TT=table;
TT.Variables=    [[{'100-250Hz'};{'100-150Hz'};{'150-200Hz'};{'200-250Hz'}] num2cell([values_spec.nl(:,1) values_spec.plusmaze(:,1) values_spec.novelty(:,1) values_spec.for(:,1) values_spec.nl(:,2) values_spec.plusmaze(:,2) values_spec.novelty(:,2) values_spec.for(:,2) values_spec.nl(:,3) values_spec.plusmaze(:,3) values_spec.novelty(:,3) values_spec.for(:,3)])];
TT.Properties.VariableNames=[{'Range'};{'HPC Baseline'};{'HPC Plusmaze'};{'HPC Novelty'};{'HPC Foraging'};{'PFC Baseline'};{'PFC Plusmaze'};{'PFC Novelty'};{'PFC Foraging'};{'PAR Baseline'};{'PAR Plusmaze'};{'PAR Novelty'};{'PAR Foraging'}];    

% % %PAR COHFOS
s=1;
w=3;
[values_spec]=getval_spectra_All(P,Q,labelconditions2,label1,s,w,win_size);
TT1=table;
TT1.Variables=    [[{'100-250Hz'};{'100-150Hz'};{'150-200Hz'};{'200-250Hz'}] num2cell([values_spec.nl(:,1) values_spec.plusmaze(:,1) values_spec.novelty(:,1) values_spec.for(:,1) values_spec.nl(:,2) values_spec.plusmaze(:,2) values_spec.novelty(:,2) values_spec.for(:,2) values_spec.nl(:,3) values_spec.plusmaze(:,3) values_spec.novelty(:,3) values_spec.for(:,3)])];
TT1.Properties.VariableNames=[{'Range'};{'HPC Baseline'};{'HPC Plusmaze'};{'HPC Novelty'};{'HPC Foraging'};{'PFC Baseline'};{'PFC Plusmaze'};{'PFC Novelty'};{'PFC Foraging'};{'PAR Baseline'};{'PAR Plusmaze'};{'PAR Novelty'};{'PAR Foraging'}];    

%HPC singles
s=2;
w=1;
[values_spec]=getval_spectra_All(P,Q,labelconditions2,label1,s,w,win_size);
TT2=table;
TT2.Variables=    [[{'100-250Hz'};{'100-150Hz'};{'150-200Hz'};{'200-250Hz'}] num2cell([values_spec.nl(:,1) values_spec.plusmaze(:,1) values_spec.novelty(:,1) values_spec.for(:,1) values_spec.nl(:,2) values_spec.plusmaze(:,2) values_spec.novelty(:,2) values_spec.for(:,2) values_spec.nl(:,3) values_spec.plusmaze(:,3) values_spec.novelty(:,3) values_spec.for(:,3)])];
TT2.Properties.VariableNames=[{'Range'};{'HPC Baseline'};{'HPC Plusmaze'};{'HPC Novelty'};{'HPC Foraging'};{'PFC Baseline'};{'PFC Plusmaze'};{'PFC Novelty'};{'PFC Foraging'};{'PAR Baseline'};{'PAR Plusmaze'};{'PAR Novelty'};{'PAR Foraging'}];    

%PAR singles
s=2;
w=3;
[values_spec]=getval_spectra_All(P,Q,labelconditions2,label1,s,w,win_size);
TT3=table;
TT3.Variables=    [[{'100-250Hz'};{'100-150Hz'};{'150-200Hz'};{'200-250Hz'}] num2cell([values_spec.nl(:,1) values_spec.plusmaze(:,1) values_spec.novelty(:,1) values_spec.for(:,1) values_spec.nl(:,2) values_spec.plusmaze(:,2) values_spec.novelty(:,2) values_spec.for(:,2) values_spec.nl(:,3) values_spec.plusmaze(:,3) values_spec.novelty(:,3) values_spec.for(:,3)])];
TT3.Properties.VariableNames=[{'Range'};{'HPC Baseline'};{'HPC Plusmaze'};{'HPC Novelty'};{'HPC Foraging'};{'PFC Baseline'};{'PFC Plusmaze'};{'PFC Novelty'};{'PFC Foraging'};{'PAR Baseline'};{'PAR Plusmaze'};{'PAR Novelty'};{'PAR Foraging'}];    


t1=repmat({'x'},[1 13]);

tab=[TT;t1;TT1;t1;TT2;t1;TT3];
%%

if win_size== 25
writetable(tab,strcat('spec_ALL_values_25_rat_', num2str(Rat),'_' ,num2str(tr(2)),'.xls'),'Sheet',1,'Range','A2:Z50')  
save('Total.mat','Total')
end

if win_size== 50
writetable(tab,strcat('spec_ALL_values_rat_', num2str(Rat),'_' ,num2str(tr(2)),'.xls'),'Sheet',1,'Range','A2:Z50')  
end
 
%% Using minimum number of ripples per type.

%HPC COHFOS
s=1;
w=1;
 n1=min([length(P.(labelconditions2{1}).(label1{w}){s}) length(P.(labelconditions2{2}).(label1{w}){s})...
        length(P.(labelconditions2{3}).(label1{w}){s}) length(P.(labelconditions2{4}).(label1{w}){s})]);
%PAR COHFOS
s=1;
w=3;
 n2=min([length(P.(labelconditions2{1}).(label1{w}){s}) length(P.(labelconditions2{2}).(label1{w}){s})...
        length(P.(labelconditions2{3}).(label1{w}){s}) length(P.(labelconditions2{4}).(label1{w}){s})]);

%HPC singles
s=2;
w=1;
 n3=min([length(P.(labelconditions2{1}).(label1{w}){s}) length(P.(labelconditions2{2}).(label1{w}){s})...
        length(P.(labelconditions2{3}).(label1{w}){s}) length(P.(labelconditions2{4}).(label1{w}){s})]);

%PAR singles
s=2;
w=3;
 n4=min([length(P.(labelconditions2{1}).(label1{w}){s}) length(P.(labelconditions2{2}).(label1{w}){s})...
        length(P.(labelconditions2{3}).(label1{w}){s}) length(P.(labelconditions2{4}).(label1{w}){s})]);
[n1 n2 n3 n4]
N=min([n1 n2 n3 n4]);

%% Find values
random_hfo=0;


win_size=50;
rand_first_run=0; %If you run for the first time.
same_nr_types=0; %Same N number across types

if same_nr_types==0
    N=[];
end
%HPC COHFOS
s=1;
w=1;
[values_spec,n1]=getval_spectra(P,Q,labelconditions2,label1,s,w,win_size,same_nr_types,N,random_hfo,rand_first_run,tr);
TT=table;
TT.Variables=    [[{'100-250Hz'};{'100-150Hz'};{'150-200Hz'};{'200-250Hz'}] num2cell([values_spec.nl(:,1) values_spec.plusmaze(:,1) values_spec.novelty(:,1) values_spec.for(:,1) values_spec.nl(:,2) values_spec.plusmaze(:,2) values_spec.novelty(:,2) values_spec.for(:,2) values_spec.nl(:,3) values_spec.plusmaze(:,3) values_spec.novelty(:,3) values_spec.for(:,3)])];
TT.Properties.VariableNames=[{'Range'};{'HPC Baseline'};{'HPC Plusmaze'};{'HPC Novelty'};{'HPC Foraging'};{'PFC Baseline'};{'PFC Plusmaze'};{'PFC Novelty'};{'PFC Foraging'};{'PAR Baseline'};{'PAR Plusmaze'};{'PAR Novelty'};{'PAR Foraging'}];    

%PAR COHFOS
s=1;
w=3;
[values_spec,n2]=getval_spectra(P,Q,labelconditions2,label1,s,w,win_size,same_nr_types,N,random_hfo,rand_first_run,tr);
TT1=table;
TT1.Variables=    [[{'100-250Hz'};{'100-150Hz'};{'150-200Hz'};{'200-250Hz'}] num2cell([values_spec.nl(:,1) values_spec.plusmaze(:,1) values_spec.novelty(:,1) values_spec.for(:,1) values_spec.nl(:,2) values_spec.plusmaze(:,2) values_spec.novelty(:,2) values_spec.for(:,2) values_spec.nl(:,3) values_spec.plusmaze(:,3) values_spec.novelty(:,3) values_spec.for(:,3)])];
TT1.Properties.VariableNames=[{'Range'};{'HPC Baseline'};{'HPC Plusmaze'};{'HPC Novelty'};{'HPC Foraging'};{'PFC Baseline'};{'PFC Plusmaze'};{'PFC Novelty'};{'PFC Foraging'};{'PAR Baseline'};{'PAR Plusmaze'};{'PAR Novelty'};{'PAR Foraging'}];    

%HPC singles
s=2;
w=1;
[values_spec,n3]=getval_spectra(P,Q,labelconditions2,label1,s,w,win_size,same_nr_types,N,random_hfo,rand_first_run,tr);
TT2=table;
TT2.Variables=    [[{'100-250Hz'};{'100-150Hz'};{'150-200Hz'};{'200-250Hz'}] num2cell([values_spec.nl(:,1) values_spec.plusmaze(:,1) values_spec.novelty(:,1) values_spec.for(:,1) values_spec.nl(:,2) values_spec.plusmaze(:,2) values_spec.novelty(:,2) values_spec.for(:,2) values_spec.nl(:,3) values_spec.plusmaze(:,3) values_spec.novelty(:,3) values_spec.for(:,3)])];
TT2.Properties.VariableNames=[{'Range'};{'HPC Baseline'};{'HPC Plusmaze'};{'HPC Novelty'};{'HPC Foraging'};{'PFC Baseline'};{'PFC Plusmaze'};{'PFC Novelty'};{'PFC Foraging'};{'PAR Baseline'};{'PAR Plusmaze'};{'PAR Novelty'};{'PAR Foraging'}];    

%PAR singles
s=2;
w=3;
[values_spec,n4]=getval_spectra(P,Q,labelconditions2,label1,s,w,win_size,same_nr_types,N,random_hfo,rand_first_run,tr);
TT3=table;
TT3.Variables=    [[{'100-250Hz'};{'100-150Hz'};{'150-200Hz'};{'200-250Hz'}] num2cell([values_spec.nl(:,1) values_spec.plusmaze(:,1) values_spec.novelty(:,1) values_spec.for(:,1) values_spec.nl(:,2) values_spec.plusmaze(:,2) values_spec.novelty(:,2) values_spec.for(:,2) values_spec.nl(:,3) values_spec.plusmaze(:,3) values_spec.novelty(:,3) values_spec.for(:,3)])];
TT3.Properties.VariableNames=[{'Range'};{'HPC Baseline'};{'HPC Plusmaze'};{'HPC Novelty'};{'HPC Foraging'};{'PFC Baseline'};{'PFC Plusmaze'};{'PFC Novelty'};{'PFC Foraging'};{'PAR Baseline'};{'PAR Plusmaze'};{'PAR Novelty'};{'PAR Foraging'}];    

[n1 n2 n3 n4]

t1=repmat({'x'},[1 13]);

tab=[TT;t1;TT1;t1;TT2;t1;TT3];
%%
if random_hfo==0

    if same_nr_types==0
        if win_size== 25
        writetable(tab,strcat('spec_values_25_rat_', num2str(Rat),'_' ,num2str(tr(2)),'.xls'),'Sheet',1,'Range','A2:Z50')  
        end

        if win_size== 50
        writetable(tab,strcat('spec_values_rat_', num2str(Rat),'_' ,num2str(tr(2)),'.xls'),'Sheet',1,'Range','A2:Z50')  
        end
    else
        if win_size== 25
        writetable(tab,strcat('spec_values_SameNR_25_rat_', num2str(Rat),'_' ,num2str(tr(2)),'.xls'),'Sheet',1,'Range','A2:Z50')  
        end

        if win_size== 50
        writetable(tab,strcat('spec_values_SameNR_rat_', num2str(Rat),'_' ,num2str(tr(2)),'.xls'),'Sheet',1,'Range','A2:Z50')  
        end
    end

else

    if same_nr_types==0
        if win_size== 25
        writetable(tab,strcat('spec_rand_values_25_rat_', num2str(Rat),'_' ,num2str(tr(2)),'.xls'),'Sheet',1,'Range','A2:Z50')  
        end

        if win_size== 50
        writetable(tab,strcat('spec_rand_values_rat_', num2str(Rat),'_' ,num2str(tr(2)),'.xls'),'Sheet',1,'Range','A2:Z50')  
        end
    else
        if win_size== 25
        writetable(tab,strcat('spec_rand_values_SameNR_25_rat_', num2str(Rat),'_' ,num2str(tr(2)),'.xls'),'Sheet',1,'Range','A2:Z50')  
        end

        if win_size== 50
        writetable(tab,strcat('spec_rand_values_SameNR_rat_', num2str(Rat),'_' ,num2str(tr(2)),'.xls'),'Sheet',1,'Range','A2:Z50')  
        end
    end

end

%% Plot spectrograms
%HPC specs
%HPC cohfos
same_nr_types=1; %Same N number across types
if same_nr_types==0
    N=[];
end
%s: 2 for single, 1 for cohfos
%w:1 for hpc centered, 3 for par centered.
s=1;
w=1;
plot_spectra(P,Q,labelconditions2,label1,s,w,same_nr_types,N)
if same_nr_types==1
printing(['Spec_HPC_cohfos_SameNR_rat' num2str(Rat) '_' num2str(tr(2))])
else
printing(['Spec_HPC_cohfos_rat' num2str(Rat) '_' num2str(tr(2))])    
end
close all


s=1;
w=3;
plot_spectra(P,Q,labelconditions2,label1,s,w,same_nr_types,N)
if same_nr_types==1
printing(['Spec_PAR_cohfos_SameNR_rat' num2str(Rat) '_' num2str(tr(2))])
else
printing(['Spec_PAR_cohfos_rat' num2str(Rat) '_' num2str(tr(2))])    
end
close all

s=2;
w=1;
plot_spectra(P,Q,labelconditions2,label1,s,w,same_nr_types,N)
if same_nr_types==1
printing(['Spec_HPC_single_SameNR_rat' num2str(Rat) '_' num2str(tr(2))])
else
printing(['Spec_HPC_single_rat' num2str(Rat) '_' num2str(tr(2))])    
end
close all


s=2;
w=3;
plot_spectra(P,Q,labelconditions2,label1,s,w,same_nr_types,N)
if same_nr_types==1
printing(['Spec_PAR_single_SameNR_rat' num2str(Rat) '_' num2str(tr(2))])
else
printing(['Spec_PAR_single_rat' num2str(Rat) '_' num2str(tr(2))])    
end
close all


%%
same_nr_types=1; %Same N number across types
if same_nr_types==0
    N=[];
end

s=1;
w=1;
plot_spec_traces(P,Q,labelconditions2,label1,s,w,same_nr_types,N)
if same_nr_types==1
printing(['SpecTraces_SameNR_HPC_cohfos_rat' num2str(Rat) '_' num2str(tr(2))])    
else
printing(['SpecTraces_HPC_cohfos_rat' num2str(Rat) '_' num2str(tr(2))])    
end
close all

s=1;
w=3;
plot_spec_traces(P,Q,labelconditions2,label1,s,w,same_nr_types,N)
if same_nr_types==1
printing(['SpecTraces_SameNR_PAR_cohfos_rat' num2str(Rat) '_' num2str(tr(2))])        
else
printing(['SpecTraces_PAR_cohfos_rat' num2str(Rat) '_' num2str(tr(2))])    
end
close all

s=2;
w=1;
plot_spec_traces(P,Q,labelconditions2,label1,s,w,same_nr_types,N)
if same_nr_types==1
printing(['SpecTraces_SameNR_HPC_single_rat' num2str(Rat) '_' num2str(tr(2))])    
else
printing(['SpecTraces_HPC_single_rat' num2str(Rat) '_' num2str(tr(2))])
end
close all

s=2;
w=3;
plot_spec_traces(P,Q,labelconditions2,label1,s,w,same_nr_types,N)
if same_nr_types==1
printing(['SpecTraces_SameNR_PAR_single_rat' num2str(Rat) '_' num2str(tr(2))])    
else
printing(['SpecTraces_PAR_single_rat' num2str(Rat) '_' num2str(tr(2))])    
end
close all




