%% GL_granger

% Main code for computation of spectral granger causality between brain areas wrt events.
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
%%
xx={'PAR'}; %Posterior Parietal cortex used to detect hfos. 

%% Get folder names
labelconditions=getfolder;
labelconditions=labelconditions.';
g=labelconditions;

multiplets=[{'singlets'} {'doublets'} {'triplets'} {'quatruplets'} {'pentuplets'} {'sextuplets'} {'septuplets'} {'octuplets'} {'nonuplets'}];
iii=1;

fn=1000; %Sampling frequency after downsampling.
%% Select conditions and sessions
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
        'Name', 'Select conditions');

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
% Load PAR signal
CORTEX=dir(strcat('*',xx{1},'*.mat'));
if isempty(CORTEX)
    g=g(~contains(g,g{k}));
    cd ..
    progress_bar(k,length(g),f)
    break
end

%PARIETAL
CORTEX=CORTEX.name;
CORTEX=load(CORTEX);
CORTEX=getfield(CORTEX,xx{1});
CORTEX=CORTEX.*(0.195); %BitVolt Open Ephys factor

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
 %Window duration in miliseconds
 ro=1200;
%Find HFOs and extraxt event-centered traces with Fieldtrip format.
[Mx_cortex,timeasleep,sig_cortex,Ex_cortex,Sx_cortex,...
  p_cortex,q_cortex,cont_cortex,sig_pq_cortex ...
  ]=gui_findripples_spec(CORTEX,states,xx,tr,PFC,HPC,fn,ro);

si=sig_cortex(~cellfun('isempty',sig_cortex));
si=[si{:}];

si_pq=sig_pq_cortex(~cellfun('isempty',sig_pq_cortex));
si_pq=[si_pq{:}];

Q_cortex=q_cortex(~cellfun('isempty',q_cortex));
Q_cortex=[Q_cortex{:}];

%Find slow and fast Hfos
[~,~,~,~,~,~,~,~,si_mixed,~]=hfo_specs(si,timeasleep,0,Rat,tr);

void_index=find(cellfun('isempty',Q_cortex));

%All par HFOS splitted in slow and fast.
Q_cortex_g1=Q_cortex(si_mixed.i1(~ismember(si_mixed.i1,void_index)));
Q_cortex_g2=Q_cortex(si_mixed.i2(~ismember(si_mixed.i2,void_index)));

%% HPC     

%Find ripples
[Mx_hpc,timeasleep,sig_hpc,Ex_hpc,Sx_hpc,...
  p_hpc,q_hpc,cont_hpc ...
]=gui_findripples_spec(HPC,states,{'HPC'},tr,PFC,CORTEX,fn,ro);

si=sig_hpc(~cellfun('isempty',sig_hpc));
si=[si{:}];

%% Coocurent events
[cohfos1,cohfos2]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),Mx_hpc,Mx_cortex,'UniformOutput',false);
%cohfos1: HPC.
%cohfos2: Cortex.

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

%%
%Cooccurrent events windows traces
p_cohfos_hpc=p_hpc(~cellfun('isempty',cohfos1));
p_cohfos_hpc=cellfun(@(equis1,equis2) equis1(equis2),p_cohfos_hpc,coh_samp_hpc,'UniformOutput',false);
p_cohfos_hpc=[p_cohfos_hpc{:}];
q_cohfos_hpc=q_hpc(~cellfun('isempty',cohfos1));
q_cohfos_hpc=cellfun(@(equis1,equis2) equis1(equis2),q_cohfos_hpc,coh_samp_hpc,'UniformOutput',false);
q_cohfos_hpc=[q_cohfos_hpc{:}];
p_cohfos_hpc=p_cohfos_hpc(~cellfun('isempty',p_cohfos_hpc));
q_cohfos_hpc=q_cohfos_hpc(~cellfun('isempty',q_cohfos_hpc));



%Single HFOs HPC
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

%Single HFOs Cortex
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
 

%% Timestamps of slow and fast HFOs
%Initialize variables
Mx_cortex_g1=Mx_cortex;
Mx_cortex_g2=Mx_cortex;
Ex_cortex_g1=Ex_cortex;
Ex_cortex_g2=Ex_cortex;
Sx_cortex_g1=Sx_cortex;
Sx_cortex_g2=Sx_cortex;

row=si_mixed.i1;
cont=0;
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

[~,cohfos2_g1]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),Mx_hpc,Mx_cortex_g1,'UniformOutput',false);
[~,cohfos2_g2]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),Mx_hpc,Mx_cortex_g2,'UniformOutput',false);

%Single HFOs Cortex
v2_g1=cellfun(@(equis1,equis2) single_hfo_get_sample(equis1,equis2),Mx_cortex_g1,cohfos2_g1,'UniformOutput',false);
%Single HFOs Cortex
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

%% Save cooccurrent, single HFOS and ripples

    %All events
     P.(strrep(labelconditions2{k},'-','_')).(label1{1})={p_cohfos_hpc;p_single_hpc;p_hpc};
     P.(strrep(labelconditions2{k},'-','_')).(label1{3})={p_cohfos_cortex;p_single_cortex;p_cortex};
    
     Q.(strrep(labelconditions2{k},'-','_')).(label1{1})={q_cohfos_hpc;q_single_hpc;q_hpc};
     Q.(strrep(labelconditions2{k},'-','_')).(label1{3})={q_cohfos_cortex;q_single_cortex;q_cortex};
     
     
     %Slow/Fast
     
     SP.(strrep(labelconditions2{k},'-','_')).(label1{3})={p_cohfos_cortex_g1;p_single_cortex_g1};
     SQ.(strrep(labelconditions2{k},'-','_')).(label1{3})={q_cohfos_cortex_g1;q_single_cortex_g1};
     
     FP.(strrep(labelconditions2{k},'-','_')).(label1{3})={p_cohfos_cortex_g2;p_single_cortex_g2};
     FQ.(strrep(labelconditions2{k},'-','_')).(label1{3})={q_cohfos_cortex_g2;q_single_cortex_g2};

progress_bar(k,length(g),f)
    cd ..    
    end
 xo
%% Granger analysis.
%P: Raw signal
%Q: Bandpassed signal.

%Labels
labelconditions3{1}='nl';
labelconditions3{2}='plusmaze';
labelconditions3{3}='novelty';
labelconditions3{4}='for';
labelconditions3=labelconditions3.';

%-------------
%Compute Granger
% s=1 cooccurrent. s=2 singles    
% w=1 HPC, w=2 PFC, w=3 PAR. Corresponds to index of 'label1'. 


s=1; %Slow Cohfos. s=1 cooccurrent. s=2 singles.
w=3; %w=3 PAR. 
[g1,g1_f,G,G_f,FB,FB1]=getval_granger(SP,SQ,labelconditions3,label1,s,w,fn);
%g1 :Parametric Granger spectrum
%g1_f :Parametric frequencies
%G :Non parametric Granger spectrum.
%G_f : Non parametric frequencies.
%FB : Frequency band granger values (Non parametric)
%FB1 : Frequency band granger values (Parametric)

%Plot Granger values across frequencies.
labelconditions3{1}='baseline';
%Parametric
granger_plot(g1,g1_f,labelconditions3,[0 300]) %All
printing(['Parametric_Slow_Cohfos' '_' num2str(tr(2))])
close all
%Non Parametric
granger_plot(G,G_f,labelconditions3,[0 300]) %All
printing(['Non_parametric_Slow_Cohfos' '_' num2str(tr(2))])
close all

%Tables with granger values for different frequency bands
TT1_20_np=table;
TT1_20_np.Variables=    [[{'PAR->PFC'};{'PFC->PAR'};{'PAR->HPC'};{'HPC->PAR'};{'PFC->HPC'};{'HPC->PFC'}] num2cell([FB{1}(:,12) FB{2}(:,12) FB{3}(:,12) FB{4}(:,12)])];
TT1_20_np.Properties.VariableNames=[{'Direction'};{labelconditions3{1}};{labelconditions3{2}};{labelconditions3{3}};{labelconditions3{4}}];    
   
TT1_20_p=table;
TT1_20_p.Variables=    [[{'PAR->PFC'};{'PFC->PAR'};{'PAR->HPC'};{'HPC->PAR'};{'PFC->HPC'};{'HPC->PFC'}] num2cell([FB1{1}(:,12) FB1{2}(:,12) FB1{3}(:,12) FB1{4}(:,12)])];
TT1_20_p.Properties.VariableNames=[{'Direction'};{labelconditions3{1}};{labelconditions3{2}};{labelconditions3{3}};{labelconditions3{4}}];    


TT1_300_np=table;
TT1_300_np.Variables=    [[{'PAR->PFC'};{'PFC->PAR'};{'PAR->HPC'};{'HPC->PAR'};{'PFC->HPC'};{'HPC->PFC'}] num2cell([FB{1}(:,13) FB{2}(:,13) FB{3}(:,13) FB{4}(:,13)])];
TT1_300_np.Properties.VariableNames=[{'Direction'};{labelconditions3{1}};{labelconditions3{2}};{labelconditions3{3}};{labelconditions3{4}}];    
   
TT1_300_p=table;
TT1_300_p.Variables=    [[{'PAR->PFC'};{'PFC->PAR'};{'PAR->HPC'};{'HPC->PAR'};{'PFC->HPC'};{'HPC->PFC'}] num2cell([FB1{1}(:,13) FB1{2}(:,13) FB1{3}(:,13) FB1{4}(:,13)])];
TT1_300_p.Properties.VariableNames=[{'Direction'};{labelconditions3{1}};{labelconditions3{2}};{labelconditions3{3}};{labelconditions3{4}}];    



%-------------
labelconditions3{1}='nl';
s=2; %Slow singles
w=3; %PAR
[g1,g1_f,G,G_f,FB,FB1,n2]=getval_granger(SP,SQ,labelconditions3,label1,s,w,fn);

labelconditions3{1}='baseline';
granger_plot(g1,g1_f,labelconditions3,[0 300]) %All
printing(['Parametric_Slow_Singles' '_' num2str(tr(2))])
close all

granger_plot(G,G_f,labelconditions3,[0 300]) %All
printing(['Non_parametric_Slow_Singles' '_' num2str(tr(2))])
close all    

%Tables with granger values for different frequency bands
TT2_20_np=table;
TT2_20_np.Variables=    [[{'PAR->PFC'};{'PFC->PAR'};{'PAR->HPC'};{'HPC->PAR'};{'PFC->HPC'};{'HPC->PFC'}] num2cell([FB{1}(:,12) FB{2}(:,12) FB{3}(:,12) FB{4}(:,12)])];
TT2_20_np.Properties.VariableNames=[{'Direction'};{labelconditions3{1}};{labelconditions3{2}};{labelconditions3{3}};{labelconditions3{4}}];    
   
TT2_20_p=table;
TT2_20_p.Variables=    [[{'PAR->PFC'};{'PFC->PAR'};{'PAR->HPC'};{'HPC->PAR'};{'PFC->HPC'};{'HPC->PFC'}] num2cell([FB1{1}(:,12) FB1{2}(:,12) FB1{3}(:,12) FB1{4}(:,12)])];
TT2_20_p.Properties.VariableNames=[{'Direction'};{labelconditions3{1}};{labelconditions3{2}};{labelconditions3{3}};{labelconditions3{4}}];    


TT2_300_np=table;
TT2_300_np.Variables=    [[{'PAR->PFC'};{'PFC->PAR'};{'PAR->HPC'};{'HPC->PAR'};{'PFC->HPC'};{'HPC->PFC'}] num2cell([FB{1}(:,13) FB{2}(:,13) FB{3}(:,13) FB{4}(:,13)])];
TT2_300_np.Properties.VariableNames=[{'Direction'};{labelconditions3{1}};{labelconditions3{2}};{labelconditions3{3}};{labelconditions3{4}}];    
   
TT2_300_p=table;
TT2_300_p.Variables=    [[{'PAR->PFC'};{'PFC->PAR'};{'PAR->HPC'};{'HPC->PAR'};{'PFC->HPC'};{'HPC->PFC'}] num2cell([FB1{1}(:,13) FB1{2}(:,13) FB1{3}(:,13) FB1{4}(:,13)])];
TT2_300_p.Properties.VariableNames=[{'Direction'};{labelconditions3{1}};{labelconditions3{2}};{labelconditions3{3}};{labelconditions3{4}}];    




%-------------    
labelconditions3{1}='nl';
s=1; %Fast Cohfos
w=3; %PAR
[g1,g1_f,G,G_f,FB,FB1,n2]=getval_granger(FP,FQ,labelconditions3,label1,s,w,fn);

labelconditions3{1}='baseline';
granger_plot(g1,g1_f,labelconditions3,[0 300]) %All
printing(['Parametric_Fast_Cohfos' '_' num2str(tr(2))])
close all

granger_plot(G,G_f,labelconditions3,[0 300]) %All
printing(['Non_parametric_Fast_Cohfos' '_' num2str(tr(2))])
close all

%Tables with granger values for different frequency bands

TT3_20_np=table;
TT3_20_np.Variables=    [[{'PAR->PFC'};{'PFC->PAR'};{'PAR->HPC'};{'HPC->PAR'};{'PFC->HPC'};{'HPC->PFC'}] num2cell([FB{1}(:,12) FB{2}(:,12) FB{3}(:,12) FB{4}(:,12)])];
TT3_20_np.Properties.VariableNames=[{'Direction'};{labelconditions3{1}};{labelconditions3{2}};{labelconditions3{3}};{labelconditions3{4}}];    
   
TT3_20_p=table;
TT3_20_p.Variables=    [[{'PAR->PFC'};{'PFC->PAR'};{'PAR->HPC'};{'HPC->PAR'};{'PFC->HPC'};{'HPC->PFC'}] num2cell([FB1{1}(:,12) FB1{2}(:,12) FB1{3}(:,12) FB1{4}(:,12)])];
TT3_20_p.Properties.VariableNames=[{'Direction'};{labelconditions3{1}};{labelconditions3{2}};{labelconditions3{3}};{labelconditions3{4}}];    


TT3_300_np=table;
TT3_300_np.Variables=    [[{'PAR->PFC'};{'PFC->PAR'};{'PAR->HPC'};{'HPC->PAR'};{'PFC->HPC'};{'HPC->PFC'}] num2cell([FB{1}(:,13) FB{2}(:,13) FB{3}(:,13) FB{4}(:,13)])];
TT3_300_np.Properties.VariableNames=[{'Direction'};{labelconditions3{1}};{labelconditions3{2}};{labelconditions3{3}};{labelconditions3{4}}];    
   
TT3_300_p=table;
TT3_300_p.Variables=    [[{'PAR->PFC'};{'PFC->PAR'};{'PAR->HPC'};{'HPC->PAR'};{'PFC->HPC'};{'HPC->PFC'}] num2cell([FB1{1}(:,13) FB1{2}(:,13) FB1{3}(:,13) FB1{4}(:,13)])];
TT3_300_p.Properties.VariableNames=[{'Direction'};{labelconditions3{1}};{labelconditions3{2}};{labelconditions3{3}};{labelconditions3{4}}];    

    
%-------------    


labelconditions3{1}='nl';    
s=2; %Fast Singles
w=3; %PAR
[g1,g1_f,G,G_f,FB,FB1,n1]=getval_granger(FP,FQ,labelconditions3,label1,s,w,fn);
labelconditions3{1}='baseline';

granger_plot(g1,g1_f,labelconditions3,[0 300]) %All
printing(['Parametric_Fast_Singles' '_' num2str(tr(2))])
close all

granger_plot(G,G_f,labelconditions3,[0 300]) %All
printing(['Non_parametric_Fast_Singles' '_' num2str(tr(2))])
close all

%Tables with granger values for different frequency bands

TT4_20_np=table;
TT4_20_np.Variables=    [[{'PAR->PFC'};{'PFC->PAR'};{'PAR->HPC'};{'HPC->PAR'};{'PFC->HPC'};{'HPC->PFC'}] num2cell([FB{1}(:,12) FB{2}(:,12) FB{3}(:,12) FB{4}(:,12)])];
TT4_20_np.Properties.VariableNames=[{'Direction'};{labelconditions3{1}};{labelconditions3{2}};{labelconditions3{3}};{labelconditions3{4}}];    
   
TT4_20_p=table;
TT4_20_p.Variables=    [[{'PAR->PFC'};{'PFC->PAR'};{'PAR->HPC'};{'HPC->PAR'};{'PFC->HPC'};{'HPC->PFC'}] num2cell([FB1{1}(:,12) FB1{2}(:,12) FB1{3}(:,12) FB1{4}(:,12)])];
TT4_20_p.Properties.VariableNames=[{'Direction'};{labelconditions3{1}};{labelconditions3{2}};{labelconditions3{3}};{labelconditions3{4}}];    


TT4_300_np=table;
TT4_300_np.Variables=    [[{'PAR->PFC'};{'PFC->PAR'};{'PAR->HPC'};{'HPC->PAR'};{'PFC->HPC'};{'HPC->PFC'}] num2cell([FB{1}(:,13) FB{2}(:,13) FB{3}(:,13) FB{4}(:,13)])];
TT4_300_np.Properties.VariableNames=[{'Direction'};{labelconditions3{1}};{labelconditions3{2}};{labelconditions3{3}};{labelconditions3{4}}];    
   
TT4_300_p=table;
TT4_300_p.Variables=    [[{'PAR->PFC'};{'PFC->PAR'};{'PAR->HPC'};{'HPC->PAR'};{'PFC->HPC'};{'HPC->PFC'}] num2cell([FB1{1}(:,13) FB1{2}(:,13) FB1{3}(:,13) FB1{4}(:,13)])];
TT4_300_p.Properties.VariableNames=[{'Direction'};{labelconditions3{1}};{labelconditions3{2}};{labelconditions3{3}};{labelconditions3{4}}];    

%-------------    


%HPC singles

labelconditions3{1}='nl';
s=2; %Slow singles
w=1;%HPC
[g1,g1_f,G,G_f,FB,FB1,n2]=getval_granger(P,Q,labelconditions3,label1,s,w,fn);

labelconditions3{1}='baseline';
granger_plot(g1,g1_f,labelconditions3,[0 300]) %All
printing(['Parametric_HPC_Singles' '_' num2str(tr(2))])
close all

granger_plot(G,G_f,labelconditions3,[0 300]) %All
printing(['Non_parametric_HPC_Singles' '_' num2str(tr(2))])
close all

%Tables with granger values for different frequency bands

TT5_20_np=table;
TT5_20_np.Variables=    [[{'PAR->PFC'};{'PFC->PAR'};{'PAR->HPC'};{'HPC->PAR'};{'PFC->HPC'};{'HPC->PFC'}] num2cell([FB{1}(:,12) FB{2}(:,12) FB{3}(:,12) FB{4}(:,12)])];
TT5_20_np.Properties.VariableNames=[{'Direction'};{labelconditions3{1}};{labelconditions3{2}};{labelconditions3{3}};{labelconditions3{4}}];    
   
TT5_20_p=table;
TT5_20_p.Variables=    [[{'PAR->PFC'};{'PFC->PAR'};{'PAR->HPC'};{'HPC->PAR'};{'PFC->HPC'};{'HPC->PFC'}] num2cell([FB1{1}(:,12) FB1{2}(:,12) FB1{3}(:,12) FB1{4}(:,12)])];
TT5_20_p.Properties.VariableNames=[{'Direction'};{labelconditions3{1}};{labelconditions3{2}};{labelconditions3{3}};{labelconditions3{4}}];    

TT5_300_np=table;
TT5_300_np.Variables=    [[{'PAR->PFC'};{'PFC->PAR'};{'PAR->HPC'};{'HPC->PAR'};{'PFC->HPC'};{'HPC->PFC'}] num2cell([FB{1}(:,13) FB{2}(:,13) FB{3}(:,13) FB{4}(:,13)])];
TT5_300_np.Properties.VariableNames=[{'Direction'};{labelconditions3{1}};{labelconditions3{2}};{labelconditions3{3}};{labelconditions3{4}}];    
   
TT5_300_p=table;
TT5_300_p.Variables=    [[{'PAR->PFC'};{'PFC->PAR'};{'PAR->HPC'};{'HPC->PAR'};{'PFC->HPC'};{'HPC->PFC'}] num2cell([FB1{1}(:,13) FB1{2}(:,13) FB1{3}(:,13) FB1{4}(:,13)])];
TT5_300_p.Properties.VariableNames=[{'Direction'};{labelconditions3{1}};{labelconditions3{2}};{labelconditions3{3}};{labelconditions3{4}}];    
%% EXPORT tables to Excel sheets

t1=repmat({'x'},[1 5]);

tab=[TT1_20_p;t1;TT2_20_p;t1;TT3_20_p;t1;TT4_20_p;t1;TT5_20_p];
writetable(tab,strcat('Granger_20_Parametric_rat_', num2str(Rat),'_' ,num2str(tr(2)),'.xls'),'Sheet',1,'Range','A2:Z50')  

tab2=[TT1_20_np;t1;TT2_20_np;t1;TT3_20_np;t1;TT4_20_np;t1;TT5_20_np];
writetable(tab2,strcat('Granger_20_NonParametric_rat_', num2str(Rat),'_' ,num2str(tr(2)),'.xls'),'Sheet',1,'Range','A2:Z50')  

tab3=[TT1_300_np;t1;TT2_300_np;t1;TT3_300_np;t1;TT4_300_np;t1;TT5_300_np];
writetable(tab3,strcat('Granger_300_NonParametric_rat_', num2str(Rat),'_' ,num2str(tr(2)),'.xls'),'Sheet',1,'Range','A2:Z50')  

tab4=[TT1_300_p;t1;TT2_300_p;t1;TT3_300_p;t1;TT4_300_p;t1;TT5_300_p];
writetable(tab4,strcat('Granger_300_Parametric_rat_', num2str(Rat),'_' ,num2str(tr(2)),'.xls'),'Sheet',1,'Range','A2:Z50')  

