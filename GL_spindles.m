%GL_spindles
%Detects spindles and computes coocurrence with hfos and ripples.
%Computes coocurrence detection for Pre and Post spindle periods.
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
%%          
%Cortical regions.
yy={'PAR'};    
xx={'PFC'};  
%Sampling freq.
fn=1000;
%% Get folder names
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
if sum(cell2mat(cellfun(@(equis1) contains(equis1,'nl'),g,'UniformOutput',false)))==1
g=g([find(cell2mat(cellfun(@(equis1) contains(equis1,labelconditions2{1}),g,'UniformOutput',false)))...
 find(cell2mat(cellfun(@(equis1) contains(equis1,labelconditions2{2}),g,'UniformOutput',false)))...
 find(cell2mat(cellfun(@(equis1) contains(equis1,labelconditions2{3}),g,'UniformOutput',false)))...
 find(cell2mat(cellfun(@(equis1) contains(equis1,labelconditions2{4}),g,'UniformOutput',false)))]);

else
    error('Name issue')
end

%Get thresholds for event detection.
tr=getfield(T,strcat('Rat',num2str(Rat)));%Thresholds 

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
CORTEX=CORTEX.*(0.195); %Open Ephys Bitvolt factor.

A = dir('*states*.mat');
A={A.name};

if  ~isempty(A)
       cellfun(@load,A);
else
      error('No Scoring found')    
end

% PFC spindles
[ripple,~,~,Mx_cortex,timeasleep,sig_cortex,Ex_cortex,Sx_cortex,...
  ~,~,~,~,~, ...
  ]=gui_findspindlesYASA(CORTEX,states,xx,multiplets,fn);


%% Cortical Ripples

CORTEX=dir(strcat('*','PAR','*.mat'));
if isempty(CORTEX)
    g=g(~contains(g,g{k}));
    cd ..
    progress_bar(k,length(g),f)
    break
end
CORTEX=CORTEX.name;
CORTEX=load(CORTEX);
CORTEX=getfield(CORTEX,'PAR');
CORTEX=CORTEX.*(0.195);

A = dir('*states*.mat');
A={A.name};
if  ~isempty(A)
       cellfun(@load,A);
else
      error('No Scoring found')    
end
%PAR HFOs detection
[ripple_sphfo,~,~,Mx_cortex_sphfo,~,sig_cortex_sphfo,Ex_cortex_sphfo,Sx_cortex_sphfo,...
  ~,~,~,~,~, ...
  ]=gui_findripples(CORTEX,states,{'PAR'},tr,multiplets,fn);

%HFOs waveforms
si=sig_cortex_sphfo(~cellfun('isempty',sig_cortex_sphfo));
si=[si{:}];
% Group in slow and fast HFOs
[~,~,~,~,~,~,~,~,si_mixed,~]=hfo_specs(si,timeasleep,0,Rat,tr);
%% Separate bimodal distribution (Average freq)  in slow and fast HFOs
%g1: slow, g2:fast
%Initialize variables
Mx_cortex_g1=Mx_cortex_sphfo;
Mx_cortex_g2=Mx_cortex_sphfo;
Ex_cortex_g1=Ex_cortex_sphfo;
Ex_cortex_g2=Ex_cortex_sphfo;
Sx_cortex_g1=Sx_cortex_sphfo;
Sx_cortex_g2=Sx_cortex_sphfo;

row=si_mixed.i1;
cont=0;
for ll=1:length(Mx_cortex_sphfo)

    if ~isempty(Mx_cortex_sphfo{ll})

        for lll=1:length(Mx_cortex_sphfo{ll})
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

%% Coocur PFC spindle and hfos
[cohfos1_g1,cohfos2_g1]=cellfun(@(equis1,equis2,equis3,equis4,equis5,equis6) co_hfo_spindle(equis1,equis2,equis3,equis4,equis5,equis6),Sx_cortex_g1,Mx_cortex_g1,Ex_cortex_g1,Sx_cortex,Mx_cortex,Ex_cortex,'UniformOutput',false);
[cohfos1_g2,cohfos2_g2]=cellfun(@(equis1,equis2,equis3,equis4,equis5,equis6) co_hfo_spindle(equis1,equis2,equis3,equis4,equis5,equis6),Sx_cortex_g2,Mx_cortex_g2,Ex_cortex_g2,Sx_cortex,Mx_cortex,Ex_cortex,'UniformOutput',false);

Cohfos1_PFC_g1{k}=([cohfos1_g1{:}]);
Cohfos2_PFC_g1{k}=([cohfos2_g1{:}]);
Cohfos1_PFC_g2{k}=([cohfos1_g2{:}]);
Cohfos2_PFC_g2{k}=([cohfos2_g2{:}]);
%%
%PRE POST ANALYSIS (PFC spindle and hfos)
%Substract/add spindle duration to compute pre/post.
[Sx_pre,Mx_pre,Ex_pre,Sx_post,Mx_post,Ex_post]=cellfun(@(equis1,equis2,equis3) pre_post_spindle(equis1,equis2,equis3) ,Sx_cortex,Mx_cortex,Ex_cortex ,'UniformOutput',false);

%Compute coocurrence for pre and post.

%PRE
[cohfos1_g1_pre,cohfos2_g1_pre]=cellfun(@(equis1,equis2,equis3,equis4,equis5,equis6) co_hfo_spindle(equis1,equis2,equis3,equis4,equis5,equis6),Sx_cortex_g1,Mx_cortex_g1,Ex_cortex_g1,Sx_pre,Mx_pre,Ex_pre,'UniformOutput',false);
[cohfos1_g2_pre,cohfos2_g2_pre]=cellfun(@(equis1,equis2,equis3,equis4,equis5,equis6) co_hfo_spindle(equis1,equis2,equis3,equis4,equis5,equis6),Sx_cortex_g2,Mx_cortex_g2,Ex_cortex_g2,Sx_pre,Mx_pre,Ex_pre,'UniformOutput',false);

Cohfos1_PFC_g1_pre{k}=([cohfos1_g1_pre{:}]);
Cohfos2_PFC_g1_pre{k}=([cohfos2_g1_pre{:}]);
Cohfos1_PFC_g2_pre{k}=([cohfos1_g2_pre{:}]);
Cohfos2_PFC_g2_pre{k}=([cohfos2_g2_pre{:}]);

%POST
[cohfos1_g1_post,cohfos2_g1_post]=cellfun(@(equis1,equis2,equis3,equis4,equis5,equis6) co_hfo_spindle(equis1,equis2,equis3,equis4,equis5,equis6),Sx_cortex_g1,Mx_cortex_g1,Ex_cortex_g1,Sx_post,Mx_post,Ex_post,'UniformOutput',false);
[cohfos1_g2_post,cohfos2_g2_post]=cellfun(@(equis1,equis2,equis3,equis4,equis5,equis6) co_hfo_spindle(equis1,equis2,equis3,equis4,equis5,equis6),Sx_cortex_g2,Mx_cortex_g2,Ex_cortex_g2,Sx_post,Mx_post,Ex_post,'UniformOutput',false);
%xo

Cohfos1_PFC_g1_post{k}=([cohfos1_g1_post{:}]);
Cohfos2_PFC_g1_post{k}=([cohfos2_g1_post{:}]);
Cohfos1_PFC_g2_post{k}=([cohfos1_g2_post{:}]);
Cohfos2_PFC_g2_post{k}=([cohfos2_g2_post{:}]);

%% Find ripples

CORTEX=dir(strcat('*','HPC','*.mat'));
if isempty(CORTEX)
    g=g(~contains(g,g{k}));
    cd ..
    progress_bar(k,length(g),f)
    break
end
CORTEX=CORTEX.name;
CORTEX=load(CORTEX);
CORTEX=getfield(CORTEX,'HPC');
CORTEX=CORTEX.*(0.195);

A = dir('*states*.mat');
A={A.name};
if  ~isempty(A)
       cellfun(@load,A);
else
      error('No Scoring found')    
end

[ripple_nhpc,~,~,Mx_nhpc,~,sig_nhpc,Ex_nhpc,Sx_nhpc,...
  ~,~,~,~,~, ...
  ]=gui_findripples(CORTEX,states,{'HPC'},tr,multiplets,fn);

%% Coocur PFC spindle and HPC ripples
[cohfos1_hpc,cohfos2_hpc]=cellfun(@(equis1,equis2,equis3,equis4,equis5,equis6) co_hfo_spindle(equis1,equis2,equis3,equis4,equis5,equis6),Sx_nhpc,Mx_nhpc,Ex_nhpc,Sx_cortex,Mx_cortex,Ex_cortex,'UniformOutput',false);


Cohfos1_PFC_hpc{k}= ([cohfos1_hpc{:}]);
Cohfos2_PFC_hpc{k}= ([cohfos2_hpc{:}]);

%PRE POST ANALYSIS
%Previously computed PRE AND POST FOR SPINDLES.

%Find coocurrent events

%PRE
[cohfos1_hpc_pre,cohfos2_hpc_pre]=cellfun(@(equis1,equis2,equis3,equis4,equis5,equis6) co_hfo_spindle(equis1,equis2,equis3,equis4,equis5,equis6),Sx_nhpc,Mx_nhpc,Ex_nhpc,Sx_pre,Mx_pre,Ex_pre,'UniformOutput',false);
Cohfos1_PFC_hpc_pre{k}= ([cohfos1_hpc_pre{:}]);
Cohfos2_PFC_hpc_pre{k}= ([cohfos2_hpc_pre{:}]);

%POST
[cohfos1_hpc_post,cohfos2_hpc_post]=cellfun(@(equis1,equis2,equis3,equis4,equis5,equis6) co_hfo_spindle(equis1,equis2,equis3,equis4,equis5,equis6),Sx_nhpc,Mx_nhpc,Ex_nhpc,Sx_post,Mx_post,Ex_post,'UniformOutput',false);
Cohfos1_PFC_hpc_post{k}= ([cohfos1_hpc_post{:}]);
Cohfos2_PFC_hpc_post{k}= ([cohfos2_hpc_post{:}]);

%% COOCUR SPINDLE-HPC MULTIPLETS
[out_pfc]=coccur_multiplets(cohfos2_hpc);
Out_PFC{k}=out_pfc;

[out_pfc_pre]=coccur_multiplets(cohfos2_hpc_pre);
Out_PFC_pre{k}=out_pfc_pre;

[out_pfc_post]=coccur_multiplets(cohfos2_hpc_post);
Out_PFC_post{k}=out_pfc_post;

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

%PAR spindles
[ripple,RipFreq,rip_duration,Mx_par,timeasleep,sig_par,Ex_par,Sx_par,...
  ripple_multiplets_par,RipFreq_multiplets_par,rip_duration_multiplets_par,sig_multiplets_par,Mx_multiplets_par...    
  ]=gui_findspindlesYASA(par,states,yy,multiplets,fn);


si=sig_par(~cellfun('isempty',sig_par));

si=[si{:}];


[x,y,z,~,~,~,l,p]=hfo_specs_spindles(si,timeasleep,fn,0);

%% Find coocur hfos and par spindles.
clear cohfos1_g1 cohfos2_g1 cohfos1_g2 cohfos2_g2 
[cohfos1_g1,cohfos2_g1]=cellfun(@(equis1,equis2,equis3,equis4,equis5,equis6) co_hfo_spindle(equis1,equis2,equis3,equis4,equis5,equis6),Sx_cortex_g1,Mx_cortex_g1,Ex_cortex_g1,Sx_par,Mx_par,Ex_par,'UniformOutput',false);
[cohfos1_g2,cohfos2_g2]=cellfun(@(equis1,equis2,equis3,equis4,equis5,equis6) co_hfo_spindle(equis1,equis2,equis3,equis4,equis5,equis6),Sx_cortex_g2,Mx_cortex_g2,Ex_cortex_g2,Sx_par,Mx_par,Ex_par,'UniformOutput',false);
%xo

Cohfos1_PAR_g1{k}=([cohfos1_g1{:}]);
Cohfos2_PAR_g1{k}=([cohfos2_g1{:}]);
Cohfos1_PAR_g2{k}=([cohfos1_g2{:}]);
Cohfos2_PAR_g2{k}=([cohfos2_g2{:}]);
%%
%PRE POST ANALYSIS
[Sx_pre,Mx_pre,Ex_pre,Sx_post,Mx_post,Ex_post]=cellfun(@(equis1,equis2,equis3) pre_post_spindle(equis1,equis2,equis3) ,Sx_par,Mx_par,Ex_par ,'UniformOutput',false);

%Find coocurrent events for PAR spindles
%PRE
[cohfos1_g1_pre,cohfos2_g1_pre]=cellfun(@(equis1,equis2,equis3,equis4,equis5,equis6) co_hfo_spindle(equis1,equis2,equis3,equis4,equis5,equis6),Sx_cortex_g1,Mx_cortex_g1,Ex_cortex_g1,Sx_pre,Mx_pre,Ex_pre,'UniformOutput',false);
[cohfos1_g2_pre,cohfos2_g2_pre]=cellfun(@(equis1,equis2,equis3,equis4,equis5,equis6) co_hfo_spindle(equis1,equis2,equis3,equis4,equis5,equis6),Sx_cortex_g2,Mx_cortex_g2,Ex_cortex_g2,Sx_pre,Mx_pre,Ex_pre,'UniformOutput',false);
%xo

Cohfos1_PAR_g1_pre{k}=([cohfos1_g1_pre{:}]);
Cohfos2_PAR_g1_pre{k}=([cohfos2_g1_pre{:}]);
Cohfos1_PAR_g2_pre{k}=([cohfos1_g2_pre{:}]);
Cohfos2_PAR_g2_pre{k}=([cohfos2_g2_pre{:}]);

%POST
[cohfos1_g1_post,cohfos2_g1_post]=cellfun(@(equis1,equis2,equis3,equis4,equis5,equis6) co_hfo_spindle(equis1,equis2,equis3,equis4,equis5,equis6),Sx_cortex_g1,Mx_cortex_g1,Ex_cortex_g1,Sx_post,Mx_post,Ex_post,'UniformOutput',false);
[cohfos1_g2_post,cohfos2_g2_post]=cellfun(@(equis1,equis2,equis3,equis4,equis5,equis6) co_hfo_spindle(equis1,equis2,equis3,equis4,equis5,equis6),Sx_cortex_g2,Mx_cortex_g2,Ex_cortex_g2,Sx_post,Mx_post,Ex_post,'UniformOutput',false);
%xo

Cohfos1_PAR_g1_post{k}=([cohfos1_g1_post{:}]);
Cohfos2_PAR_g1_post{k}=([cohfos2_g1_post{:}]);
Cohfos1_PAR_g2_post{k}=([cohfos1_g2_post{:}]);
Cohfos2_PAR_g2_post{k}=([cohfos2_g2_post{:}]);


%% Coocur PAR spindle and HPC ripples
clear cohfos1_hpc cohfos2_hpc 
[cohfos1_hpc,cohfos2_hpc]=cellfun(@(equis1,equis2,equis3,equis4,equis5,equis6) co_hfo_spindle(equis1,equis2,equis3,equis4,equis5,equis6),Sx_nhpc,Mx_nhpc,Ex_nhpc,Sx_par,Mx_par,Ex_par,'UniformOutput',false);


Cohfos1_PAR_hpc{k}= ([cohfos1_hpc{:}]);
Cohfos2_PAR_hpc{k}= ([cohfos2_hpc{:}]);

%PRE POST ANALYSIS
%PRE AND POST SPINDLE COMPUTED ABOVE.

%PRE
[cohfos1_hpc_pre,cohfos2_hpc_pre]=cellfun(@(equis1,equis2,equis3,equis4,equis5,equis6) co_hfo_spindle(equis1,equis2,equis3,equis4,equis5,equis6),Sx_nhpc,Mx_nhpc,Ex_nhpc,Sx_pre,Mx_pre,Ex_pre,'UniformOutput',false);
Cohfos1_PAR_hpc_pre{k}= ([cohfos1_hpc_pre{:}]);
Cohfos2_PAR_hpc_pre{k}= ([cohfos2_hpc_pre{:}]);

%POST
[cohfos1_hpc_post,cohfos2_hpc_post]=cellfun(@(equis1,equis2,equis3,equis4,equis5,equis6) co_hfo_spindle(equis1,equis2,equis3,equis4,equis5,equis6),Sx_nhpc,Mx_nhpc,Ex_nhpc,Sx_post,Mx_post,Ex_post,'UniformOutput',false);
Cohfos1_PAR_hpc_post{k}= ([cohfos1_hpc_post{:}]);
Cohfos2_PAR_hpc_post{k}= ([cohfos2_hpc_post{:}]);


%% COOCUR PAR SPINDLE-HPC MULTIPLETS
[out_par]=coccur_multiplets(cohfos2_hpc);
Out_PAR{k}=out_par;

[out_par_pre]=coccur_multiplets(cohfos2_hpc_pre);
Out_PAR_pre{k}=out_par_pre;

[out_par_post]=coccur_multiplets(cohfos2_hpc_post);
Out_PAR_post{k}=out_par_post;

%% Co-occurrent spindles
[cohfos1,cohfos2]=cellfun(@(equis1,equis2,equis3,equis4,equis5,equis6) co_hfo_spindle(equis1,equis2,equis3,equis4,equis5,equis6),Sx_par,Mx_par,Ex_par,Sx_cortex,Mx_cortex,Ex_cortex,'UniformOutput',false);
%Remove repeated values
[cohfos1,cohfos2]=cellfun(@(equis1,equis2) coocur_repeat(equis1,equis2), cohfos1,cohfos2,'UniformOutput',false);

%cohfos1: PAR.
%cohfos2: PFC.

%% Co-occurrent spindles
cohf_mx_par=Mx_par(~cellfun('isempty',cohfos1));%Peak values cells where co-occurring events were found.
cohf_sx_par=Sx_par(~cellfun('isempty',cohfos1));%Peak values cells where co-occurring events were found.
cohf_ex_par=Ex_par(~cellfun('isempty',cohfos1));%Peak values cells where co-occurring events were found.

Cohfos1=cohfos1(~cellfun('isempty',cohfos1));

%Locate sample per cohfos
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


%Single PAR 
v2=cellfun(@(equis1,equis2) single_hfo_get_sample(equis1,equis2),Mx_par,cohfos1,'UniformOutput',false);

Sig_par_single=cellfun(@(equis1,equis2) equis1(equis2),sig_par,v2,'UniformOutput',false);
Sig_par_single=[Sig_par_single{:}];


[single_mx_par_val,single_sx_par_val]=cellfun(@(equis1,equis2,equis3) single_hfos_mx(equis1,equis2,equis3),cohfos1,Mx_par,Sx_par,'UniformOutput',false);
single_mx_par_val=[single_mx_par_val{:}];
single_sx_par_val=[single_sx_par_val{:}];


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


%Single HFOs Cortex
v2=cellfun(@(equis1,equis2) single_hfo_get_sample(equis1,equis2),Mx_cortex,cohfos2,'UniformOutput',false);

Sig_cortex_single=cellfun(@(equis1,equis2) equis1(equis2),sig_cortex,v2,'UniformOutput',false);
Sig_cortex_single=[Sig_cortex_single{:}];

[single_mx_cortex_val,single_sx_cortex_val]=cellfun(@(equis1,equis2,equis3) single_hfos_mx(equis1,equis2,equis3),cohfos2,Mx_cortex,Sx_cortex,'UniformOutput',false);
single_mx_cortex_val=[single_mx_cortex_val{:}];
single_sx_cortex_val=[single_sx_cortex_val{:}];

progress_bar(k,length(g),f)
    cd ..    
    end
  xo
%%

save('Cohfos1_PAR_g1_pre','Cohfos1_PAR_g2_pre','Cohfos1_PFC_g1_pre','Cohfos1_PFC_g2_pre',...
    'Cohfos1_PAR_g1_post','Cohfos1_PAR_g2_post','Cohfos1_PFC_g1_post','Cohfos1_PFC_g2_post',...
    'Cohfos2_PAR_hpc_pre','Cohfos2_PAR_hpc_post','Cohfos2_PFC_hpc_pre','Cohfos2_PFC_hpc_post',...
    'Out_PAR_pre','Out_PFC_pre','Out_PAR_post','Out_PFC_post');

  %% HFOs-spindles PRE and POST
    TT=table;
    TT.Variables= [cellfun('length',Cohfos1_PAR_g1_pre); cellfun('length',Cohfos1_PAR_g2_pre);...
       cellfun('length',Cohfos1_PFC_g1_pre); cellfun('length',Cohfos1_PFC_g2_pre); 
    ] ;
    TT.Properties.VariableNames=[cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false).'].';
    writetable(TT,strcat('ripple_spindles_pre','.xls'),'Sheet',1,'Range','A2:L10')    

    TT=table;
    TT.Variables= [cellfun('length',Cohfos1_PAR_g1_post); cellfun('length',Cohfos1_PAR_g2_post);...
       cellfun('length',Cohfos1_PFC_g1_post); cellfun('length',Cohfos1_PFC_g2_post); 
    ] ;
    TT.Properties.VariableNames=[cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false).'].';
    writetable(TT,strcat('ripple_spindles_post','.xls'),'Sheet',1,'Range','A2:L10')    

%% HPC ripples-spindles PRE POST
    TT=table;
    TT.Variables= [cellfun('length',Cohfos2_PAR_hpc_pre); cell2mat(cellfun(@(equis) length(unique(equis)),Cohfos2_PAR_hpc_pre,'UniformOutput',false));...
       cellfun('length',Cohfos2_PFC_hpc_pre); cell2mat(cellfun(@(equis) length(unique(equis)),Cohfos2_PFC_hpc_pre,'UniformOutput',false))]; 
    
    TT.Properties.VariableNames=[cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false).'].';
                writetable(TT,strcat('hpc_ripple_spindles_pre','.xls'),'Sheet',1,'Range','A2:L50')    

                    TT=table;
    TT.Variables= [cellfun('length',Cohfos2_PAR_hpc_post); cell2mat(cellfun(@(equis) length(unique(equis)),Cohfos2_PAR_hpc_post,'UniformOutput',false));...
       cellfun('length',Cohfos2_PFC_hpc_post); cell2mat(cellfun(@(equis) length(unique(equis)),Cohfos2_PFC_hpc_post,'UniformOutput',false))]; 
    
    TT.Properties.VariableNames=[cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false).'].';
                writetable(TT,strcat('hpc_ripple_spindles_post','.xls'),'Sheet',1,'Range','A2:L50')    

%% HPC ripples-spindles-multiplets PRE POST
    TT=table;
    TT.Variables=([Out_PAR_pre{1}(:,2) Out_PAR_pre{2}(:,2) Out_PAR_pre{3}(:,2) Out_PAR_pre{4}(:,2); [NaN NaN NaN NaN] ;Out_PFC_pre{1}(:,2) Out_PFC_pre{2}(:,2) Out_PFC_pre{3}(:,2) Out_PFC_pre{4}(:,2)] );
    TT.Properties.VariableNames=[cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false).'].';
                 writetable(TT,strcat('hpc_ripple_spindles_multiplets_pre','.xls'),'Sheet',1,'Range','A2:L50')    

    TT=table;
    TT.Variables=([Out_PAR_post{1}(:,2) Out_PAR_post{2}(:,2) Out_PAR_post{3}(:,2) Out_PAR_post{4}(:,2); [NaN NaN NaN NaN] ;Out_PFC_post{1}(:,2) Out_PFC_post{2}(:,2) Out_PFC_post{3}(:,2) Out_PFC_post{4}(:,2)] );
    TT.Properties.VariableNames=[cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false).'].';
                 writetable(TT,strcat('hpc_ripple_spindles_multiplets_post','.xls'),'Sheet',1,'Range','A2:L50')    

  %% HFOs-spindles
    TT=table;
    TT.Variables= [cellfun('length',Cohfos1_PAR_g1); cellfun('length',Cohfos1_PAR_g2);...
       cellfun('length',Cohfos1_PFC_g1); cellfun('length',Cohfos1_PFC_g2); 
    ] ;
    TT.Properties.VariableNames=[cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false).'].';
            writetable(TT,strcat('ripple_spindles','.xls'),'Sheet',1,'Range','A2:L10')    
%% HPC ripples-spindles
    TT=table;
    TT.Variables= [cellfun('length',Cohfos2_PAR_hpc); cell2mat(cellfun(@(equis) length(unique(equis)),Cohfos2_PAR_hpc,'UniformOutput',false));...
       cellfun('length',Cohfos2_PFC_hpc); cell2mat(cellfun(@(equis) length(unique(equis)),Cohfos2_PFC_hpc,'UniformOutput',false))]; 
    
    TT.Properties.VariableNames=[cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false).'].';
                writetable(TT,strcat('hpc_ripple_spindles','.xls'),'Sheet',1,'Range','A2:L50')    
%% HPC ripples-spindles-multiplets
    TT=table;
TT.Variables=([Out_PAR{1}(:,2) Out_PAR{2}(:,2) Out_PAR{3}(:,2) Out_PAR{4}(:,2); [NaN NaN NaN NaN] ;Out_PFC{1}(:,2) Out_PFC{2}(:,2) Out_PFC{3}(:,2) Out_PFC{4}(:,2)] );
    TT.Properties.VariableNames=[cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false).'].';
                writetable(TT,strcat('hpc_ripple_spindles_multiplets','.xls'),'Sheet',1,'Range','A2:L50')    

