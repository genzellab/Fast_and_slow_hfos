%GL_spindles_counts.
%Detects spindles and computes coocurrence with hfos and ripples.
% Requires 'load_me_first.mat' loaded first. 

%% Find location
close all
dname=uigetdir([],'Select folder with Matlab data containing all rats.');
cd(dname)
%cd('/home/adrian/Documents/Plusmaze_downsampled')

%%
%Select rat number
opts.Resize = 'on';
opts.WindowStyle = 'modal';
opts.Interpreter = 'tex';
prompt=strcat('\bf Select a rat#. Options:','{ }',num2str(rats));
answer = inputdlg(prompt,'Input',[2 30],{''},opts);
Rat=str2num(answer{1});
cd(num2str(Rat))
%%
%Cortical regions.
yy={'PAR'};    
xx={'PFC'};  
%Sampling freq.
fn=1000;

%%
labelconditions=getfolder;
labelconditions=labelconditions.';
g=labelconditions;

multiplets=[{'singlets'} {'doublets'} {'triplets'} {'quatruplets'} {'pentuplets'} {'sextuplets'} {'septuplets'} {'octuplets'} {'nonuplets'}];
iii=1;

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
%xo
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
CORTEX=dir(strcat('*',xx{1},'*.mat'));
if isempty(CORTEX)
    g=g(~contains(g,g{k}));
    cd ..
    progress_bar(k,length(g),f)
    break
end
CORTEX=CORTEX.name;
CORTEX=load(CORTEX);
CORTEX=getfield(CORTEX,xx{1});
CORTEX=CORTEX.*(0.195);

A = dir('*states*.mat');
A={A.name};

if  ~isempty(A)
       cellfun(@load,A);
else
      error('No Scoring found')    
end

%Find PFC spindles
[ripple,RipFreq,rip_duration,Mx_cortex,timeasleep,sig_cortex,Ex_cortex,Sx_cortex,...
  ripple_multiplets_cortex,RipFreq_multiplets_cortex,rip_duration_multiplets_cortex,sig_multiplets_cortex,~, ...
  ]=gui_findspindlesYASA(CORTEX,states,xx,multiplets,fn);

si=sig_cortex(~cellfun('isempty',sig_cortex));
si=[si{:}];

All_Par.( strrep(g{k},'-','_'))=si;

print_hist=0; %Print histogram? yes=1, no=0.
[x,y,z,~,~,~,l,p]=hfo_specs_spindles(si,timeasleep,fn,print_hist);
if print_hist==1
cd ..
printing(['Histograms_Spindles_' xx{1} '_Count_' g{k}]);
close all
cd(g{k})
end

fi_cortex(k)=x;
fa_cortex(k)=y;
amp_cortex(k)=z;
auc_cortex(k)=l;
p2p_cortex(k)=p;

%% PFC spindles
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
    
%% PARIETAL SPINDLES. 

par=dir(strcat('*',yy{1},'*.mat'));
if isempty(par)
    g=g(~contains(g,g{k}));
    cd ..
    progress_bar(k,length(g),f)
    break
end

par=par.name;
par=load(par);
par=getfield(par,yy{1});
par=par.*(0.195);

[ripple,RipFreq,rip_duration,Mx_par,timeasleep,sig_par,Ex_par,Sx_par,...
  ripple_multiplets_par,RipFreq_multiplets_par,rip_duration_multiplets_par,sig_multiplets_par,Mx_multiplets_par...    
  ]=gui_findspindlesYASA(par,states,yy,multiplets,fn);

si=sig_par(~cellfun('isempty',sig_par));
si=[si{:}];

All_par.( strrep(g{k},'-','_'))=si;

print_hist=0;
[x,y,z,~,~,~,l,p]=hfo_specs_spindles(si,timeasleep,fn,print_hist);
if print_hist==1
    cd ..
    printing(['Histograms_Spindles_' yy{1} '_Count_' g{k}]);
    close all
    cd(g{k})
end
%Main features of spindle
fi_par(k)=x;
fa_par(k)=y;
amp_par(k)=z;
auc_par(k)=l;
p2p_par(k)=p;

%% PAR spindles
hfos_par(k)=ripple;
hfos_par_rate(k)=RipFreq;
hfos_par_duration(k)=rip_duration;

%Multiplets    
for ll=1:length(multiplets)
   eval(['hfos_par_' multiplets{ll} '(k)=ripple_multiplets_par.' multiplets{ll} ';']) 
   eval(['hfos_par_rate_' multiplets{ll} '(k)=RipFreq_multiplets_par.' multiplets{ll} ';']) 
   eval(['hfos_par_duration_' multiplets{ll} '(k)=rip_duration_multiplets_par.' multiplets{ll} ';'])    
end
% xo
%% Coocurent spindles
[cohfos1,cohfos2]=cellfun(@(equis1,equis2,equis3,equis4,equis5,equis6) co_hfo_spindle(equis1,equis2,equis3,equis4,equis5,equis6),Sx_par,Mx_par,Ex_par,Sx_cortex,Mx_cortex,Ex_cortex,'UniformOutput',false);
%Remove repeated values
[cohfos1,cohfos2]=cellfun(@(equis1,equis2) coocur_repeat(equis1,equis2), cohfos1,cohfos2,'UniformOutput',false);


%cohfos1: par.
%cohfos2: Cortex.
%Common values:
cohfos_count(k)=sum(cellfun('length',cohfos1));
cohfos_rate(k)=sum(cellfun('length',cohfos1))/(timeasleep*(60));

%Multiplet cohfos
for ll=1:length(multiplets)
[cohfos1_multiplets.(multiplets{ll}),cohfos2_multiplets.(multiplets{ll})]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),Mx_multiplets_par.(multiplets{ll}).',Mx_cortex,'UniformOutput',false);
cohfos_count_multiplets.(multiplets{ll})(k)=sum(cellfun('length',cohfos1_multiplets.(multiplets{ll})));
cohfos_rate_multiplets.(multiplets{ll})(k)=sum(cellfun('length',cohfos1_multiplets.(multiplets{ll})))/(timeasleep*(60));
end
%%

%par COHFOS
cohf_mx_par=Mx_par(~cellfun('isempty',cohfos1));%Peak values cells where co-occurrent events were found.
cohf_sx_par=Sx_par(~cellfun('isempty',cohfos1));%Peak values cells where co-occurrent events were found.
cohf_ex_par=Ex_par(~cellfun('isempty',cohfos1));%Peak values cells where co-occurrent events were found.

Cohfos1=cohfos1(~cellfun('isempty',cohfos1));

%Locate sample per co-occurrent events
coh_samp_par= cellfun(@(equis1,equis2) co_hfo_get_sample(equis1,equis2),cohf_mx_par,Cohfos1,'UniformOutput',false);

cohf_sx_par_val=cellfun(@(equis1,equis2) equis1(equis2),cohf_sx_par,coh_samp_par,'UniformOutput',false);
cohf_sx_par_val=[cohf_sx_par_val{:}];

cohf_mx_par_val=cellfun(@(equis1,equis2) equis1(equis2),cohf_mx_par,coh_samp_par,'UniformOutput',false);
cohf_mx_par_val=[cohf_mx_par_val{:}];

cohf_ex_par_val=cellfun(@(equis1,equis2) equis1(equis2),cohf_ex_par,coh_samp_par,'UniformOutput',false);
cohf_ex_par_val=[cohf_ex_par_val{:}];

cohf_par_dura=cohf_ex_par_val-cohf_sx_par_val;
cohf_par_dura=median(cohf_par_dura);
Cohf_par_dura(k)=cohf_par_dura;

Sig_par=sig_par(~cellfun('isempty',cohfos1));
Sig_par=cellfun(@(equis1,equis2) equis1(equis2),Sig_par,coh_samp_par,'UniformOutput',false);
Sig_par=[Sig_par{:}];

[x,y,z,w,h,q,l,p]=hfo_specs_spindles(Sig_par,timeasleep,fn,0);
fi_cohfo_par(k)=x;
fa_cohfo_par(k)=y;
amp_cohfo_par(k)=z;
count_cohfo_par(k)=w;
rate_cohfo_par(k)=h;
dura_cohfo_par(k)=q;
auc_cohfo_par(k)=l;
p2p_cohfo_par(k)=p;

%Single events
v2=cellfun(@(equis1,equis2) single_hfo_get_sample(equis1,equis2),Mx_par,cohfos1,'UniformOutput',false);

Sig_par_single=cellfun(@(equis1,equis2) equis1(equis2),sig_par,v2,'UniformOutput',false);
Sig_par_single=[Sig_par_single{:}];


[single_mx_par_val,single_sx_par_val]=cellfun(@(equis1,equis2,equis3) single_hfos_mx(equis1,equis2,equis3),cohfos1,Mx_par,Sx_par,'UniformOutput',false);
single_mx_par_val=[single_mx_par_val{:}];
single_sx_par_val=[single_sx_par_val{:}];

[x,y,z,w,h,q,l,p]=hfo_specs_spindles(Sig_par_single,timeasleep,fn,0);
fi_single_par(k)=x;
fa_single_par(k)=y;
amp_single_par(k)=z;
count_single_par(k)=w;
rate_single_par(k)=h;
dura_single_par(k)=q;
auc_single_par(k)=l;
p2p_single_par(k)=p;


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
cohf_cortex_dura=cohf_ex_cortex_val-cohf_sx_cortex_val;
cohf_cortex_dura=median(cohf_cortex_dura);
Cohf_cortex_dura(k)=cohf_cortex_dura;

Sig_cortex=sig_cortex(~cellfun('isempty',cohfos2));
Sig_cortex=cellfun(@(equis1,equis2) equis1(equis2),Sig_cortex,coh_samp_cortex,'UniformOutput',false);
Sig_cortex=[Sig_cortex{:}];

[x,y,z,w,h,q,l,p]=hfo_specs_spindles(Sig_cortex,timeasleep,fn,0);
fi_cohfo_cortex(k)=x;
fa_cohfo_cortex(k)=y;
amp_cohfo_cortex(k)=z;
count_cohfo_cortex(k)=w;
rate_cohfo_cortex(k)=h;
dura_cohfo_cortex(k)=q;
auc_cohfo_cortex(k)=l;
p2p_cohfo_cortex(k)=p;

%Single events Cortex
v2=cellfun(@(equis1,equis2) single_hfo_get_sample(equis1,equis2),Mx_cortex,cohfos2,'UniformOutput',false);

Sig_cortex_single=cellfun(@(equis1,equis2) equis1(equis2),sig_cortex,v2,'UniformOutput',false);
Sig_cortex_single=[Sig_cortex_single{:}];

[single_mx_cortex_val,single_sx_cortex_val]=cellfun(@(equis1,equis2,equis3) single_hfos_mx(equis1,equis2,equis3),cohfos2,Mx_cortex,Sx_cortex,'UniformOutput',false);
single_mx_cortex_val=[single_mx_cortex_val{:}];
single_sx_cortex_val=[single_sx_cortex_val{:}];


[x,y,z,w,h,q,l,p]=hfo_specs_spindles(Sig_cortex_single,timeasleep,fn,0);
fi_single_cortex(k)=x;
fa_single_cortex(k)=y;
amp_single_cortex(k)=z;
count_single_cortex(k)=w;
rate_single_cortex(k)=h;
dura_single_cortex(k)=q;
auc_single_cortex(k)=l;
p2p_single_cortex(k)=p;

progress_bar(k,length(g),f)
    cd ..    
    end
  xo
%Plot histogram of main features for PFC and PAR spindles

A_cell = struct2cell(All_Par);
All_Par_26=[A_cell{:}];
hfo_specs_spindles(All_Par_26,timeasleep,fn,1)
printing(['HistogramSpindles_' num2str(Rat) '_All_PFC_count'])
close all


A_cell = struct2cell(All_par);
All_Par_26=[A_cell{:}];
hfo_specs_spindles(All_Par_26,timeasleep,fn,1)
printing(['HistogramSpindles_' num2str(Rat) '_All_PAR_count'])
close all

%% Counts of spindles, cooccurrent spindles and single spindles.
    TT=table;
    TT.Variables=    [[{['Count_' xx{1} '_total']};{['Count_' yy{1} '_total']};{['Count_coocur_' xx{1} '_' yy{1}]};{['Count_single_' xx{1}]};{['Count_single_' yy{1}]}] num2cell([...
        hfos_cortex;hfos_par;count_cohfo_par;count_single_cortex;count_single_par])];

    TT.Properties.VariableNames=['Metric';cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false).'].';

            writetable(TT,strcat(xx{1},'_',yy{1},'_spindles','.xls'),'Sheet',1,'Range','A2:L10')    
%% PFC spindles features
    %All detections cortex
    TT=table;
    TT.Variables=    [[{'Count'};{'Rate'};{'Duration'};{'Average Frequency'};{'Instantaneous Frequency'};{'Amplitude'}] num2cell([hfos_cortex;hfos_cortex_rate;hfos_cortex_duration;fa_cortex;fi_cortex; amp_cortex])];    
    TT.Properties.VariableNames=['Metric' cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)].';
    
            writetable(TT,strcat(xx{1},'_Features_spindles','.xls'),'Sheet',1,'Range','A2:L10')    

%%
%Cortex coocur features
    TT=table;
    TT.Variables=    [[{'Count'};{'Rate'};{'Duration'};{'Average Frequency'};{'Instantaneous Frequency'};{'Amplitude'}] num2cell([count_cohfo_cortex;rate_cohfo_cortex;dura_cohfo_cortex;fa_cohfo_cortex;fi_cohfo_cortex; amp_cohfo_cortex])];
    
    TT.Properties.VariableNames=['Metric' cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)].';
    
            writetable(TT,strcat(xx{1},'_spindles','_cohfos','.xls'),'Sheet',1,'Range','A2:L10')    

%Cortex singles features
    TT=table;
    TT.Variables=    [[{'Count'};{'Rate'};{'Duration'};{'Average Frequency'};{'Instantaneous Frequency'};{'Amplitude'}] num2cell([count_single_cortex;rate_single_cortex;dura_single_cortex;fa_single_cortex;fi_single_cortex; amp_single_cortex])];
    TT.Properties.VariableNames=['Metric' cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)].';
    writetable(TT,strcat(xx{1},'_spindles','_singles','.xls'),'Sheet',1,'Range','A2:L10')    

  %%  PAR spindles features      
%All detections
    TT=table;
    TT.Variables=    [[{'Count'};{'Rate'};{'Duration'};{'Average Frequency'};{'Instantaneous Frequency'};{'Amplitude'}] num2cell([hfos_par;hfos_par_rate;hfos_par_duration;fa_par;fi_par;amp_par])];        
    TT.Properties.VariableNames=['Metric' cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)].';
    
            writetable(TT,strcat(yy{1},'_Features_spindles','.xls'),'Sheet',1,'Range','A2:L10')    


%%
%PAR spindles coocur
    TT=table;
    TT.Variables=    [[{'Count'};{'Rate'};{'Duration'};{'Average Frequency'};{'Instantaneous Frequency'};{'Amplitude'}] num2cell([count_cohfo_par;rate_cohfo_par;dura_cohfo_par;fa_cohfo_par;fi_cohfo_par; amp_cohfo_par])];
    
    TT.Properties.VariableNames=['Metric' cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)].';
    
            writetable(TT,strcat(yy{1},'_spindles','_cohfos_','.xls'),'Sheet',1,'Range','A2:L10')    

%PAR spindles singles
    TT=table;
    TT.Variables=    [[{'Count'};{'Rate'};{'Duration'};{'Average Frequency'};{'Instantaneous Frequency'};{'Amplitude'}] num2cell([count_single_par;rate_single_par;dura_single_par;fa_single_par;fi_single_par; amp_single_par])];
    
    TT.Properties.VariableNames=['Metric' cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)].';
    
            writetable(TT,strcat(yy{1},'_spindles','_singles','.xls'),'Sheet',1,'Range','A2:L10')    

