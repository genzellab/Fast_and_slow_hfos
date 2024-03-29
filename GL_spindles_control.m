%GL_spindles_control
%Detects spindles and computes coocurrence with hfos and ripples.
%Controls by changing timestamps randomly to generate a null-distribution.
% Requires 'load_me_first.mat' loaded first. 

%% Find location
% close all
% dname=uigetdir([],'Select folder with Matlab data containing all rats.');
% cd(dname)
cd('/home/adrian/Documents/Plusmaze_downsampled')

%%
%Select rat ID
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

%%  Select conditions and sessions
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
CORTEX=CORTEX.*(0.195);

A = dir('*states*.mat');
A={A.name};

if  ~isempty(A)
       cellfun(@load,A);
else
      error('No Scoring found')    
end

%Find PFC Spindles
clear Ex_cortex Sx_cortex Mx_cortex
[ripple,~,~,Mx_cortex,timeasleep,sig_cortex,Ex_cortex,Sx_cortex,...
  ~,~,~,~,~, ...
  ]=gui_findspindlesYASA(CORTEX,states,xx,multiplets,fn);

si=sig_cortex(~cellfun('isempty',sig_cortex));
si=[si{:}];

%PRE POST ANALYSIS
% [Sx_pre,Mx_pre,Ex_pre,Sx_post,Mx_post,Ex_post]=cellfun(@(equis1,equis2,equis3) pre_post_spindle(equis1,equis2,equis3) ,Sx_cortex,Mx_cortex,Ex_cortex ,'UniformOutput',false);

%%
%Convert signal to 1 sec epochs.
        e_t=1;
        e_samples=e_t*(fn); %fs=1kHz
        ch=length(CORTEX);
        nc=floor(ch/e_samples); %Number of epochsw
        NC=[];
        for kk=1:nc
          NC(:,kk)= CORTEX(1+e_samples*(kk-1):e_samples*kk);
        end
        %IF NREM vec_bin=1.
        vec_bin=states;
        vec_bin(vec_bin~=3)=0;
        vec_bin(vec_bin==3)=1;
        %Cluster one values:
        v2=ConsecutiveOnes(vec_bin);
        v_index=find(v2~=0);
        v_values=v2(v2~=0);
        clear v
    for epoch_count=1:length(v_index)
    v{epoch_count,1}=reshape(NC(:, v_index(epoch_count):v_index(epoch_count)+(v_values(1,epoch_count)-1)), [], 1);
    end    
    ti=cellfun(@(equis) reshape(linspace(0, length(equis)-1,length(equis))*(1/fn),[],1) ,v,'UniformOutput',false);

%%
clear dur_cortex
%Find duration of spindles
for dd=1:length(Mx_cortex)
        dur_cortex{dd}=Ex_cortex{dd}-Sx_cortex{dd};
    if isempty(Sx_cortex{dd})
        Sx_cortex{dd}=nan;
    end

end
dur_cortex=dur_cortex.';
    clear Sr_cortex Er_cortex
% Shuffling timestamps randomly 1000 times    
    rng('default') %Initialize seed for repeatability
    parfor r=1:1000
        ti_rand=cellfun(@(equis) equis(randperm(size(equis, 1))),ti,'UniformOutput',false);
        Sr_cortex{r}=cellfun(@(equis1,equis2,equis3) equis3(findclosest(equis1,equis2)).', ti,Sx_cortex,ti_rand,'UniformOutput',false );
        Er_cortex{r}=cellfun(@(equis1,equis2,equis3,equis4) equis3(findclosest(equis1,equis2)).'+equis4, ti,Sx_cortex,ti_rand,dur_cortex,'UniformOutput',false );
    end

%Change Nans to []
for dd=1:length(Mx_cortex)
    if isnan(Sx_cortex{dd})
        Sx_cortex{dd}=[];
    end

end

%% High Frequency oscillations

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
% Detect HFOs
[ripple_sphfo,~,~,Mx_cortex_sphfo,~,sig_cortex_sphfo,Ex_cortex_sphfo,Sx_cortex_sphfo,...
  ~,~,~,~,~, ...
  ]=gui_findripples(CORTEX,states,{'PAR'},tr,multiplets,fn);


si=sig_cortex_sphfo(~cellfun('isempty',sig_cortex_sphfo));
si=[si{:}];

%Group events in slow and fast HFOs.
[~,~,~,~,~,~,~,~,si_mixed,~]=hfo_specs(si,timeasleep,0,Rat,tr);
%% Determining slow and fast HFOs timestamps.
%Initializing variables
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
%% Coocur PFC spindle and shuffled hfos.

parfor r=1:1000
[cohfos1_rand_g1,cohfos2_rand_g1]=cellfun(@(equis1,equis2,equis3,equis4,equis5,equis6) co_hfo_spindle(equis1,equis2,equis3,equis4,equis5,equis6),Sx_cortex_g1,Mx_cortex_g1,Ex_cortex_g1,Sr_cortex{r},Sr_cortex{r},Er_cortex{r},'UniformOutput',false);
[cohfos1_rand_g2,cohfos2_rand_g2]=cellfun(@(equis1,equis2,equis3,equis4,equis5,equis6) co_hfo_spindle(equis1,equis2,equis3,equis4,equis5,equis6),Sx_cortex_g2,Mx_cortex_g2,Ex_cortex_g2,Sr_cortex{r},Sr_cortex{r},Er_cortex{r},'UniformOutput',false);


Cohfos1_PFC_all_g1(k,r)=sum(cellfun('length',cohfos1_rand_g1));
Cohfos1_PFC_unique_g1(k,r)=sum(cellfun(@(equis) length(unique(equis)), cohfos1_rand_g1));
Cohfos2_PFC_all_g1(k,r)=sum(cellfun('length',cohfos2_rand_g1));
Cohfos2_PFC_unique_g1(k,r)=sum(cellfun(@(equis) length(unique(equis)), cohfos2_rand_g1));

Cohfos1_PFC_all_g2(k,r)=sum(cellfun('length',cohfos1_rand_g2));
Cohfos1_PFC_unique_g2(k,r)=sum(cellfun(@(equis) length(unique(equis)), cohfos1_rand_g2));
Cohfos2_PFC_all_g2(k,r)=sum(cellfun('length',cohfos2_rand_g2));
Cohfos2_PFC_unique_g2(k,r)=sum(cellfun(@(equis) length(unique(equis)), cohfos2_rand_g2));
end

%% HPC ripples

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
% Find ripples
[ripple_nhpc,~,~,Mx_nhpc,~,sig_nhpc,Ex_nhpc,Sx_nhpc,...
  ~,~,~,~,~, ...
  ]=gui_findripples(CORTEX,states,{'HPC'},tr,multiplets,fn);

%% Coocur PFC spindle and HPC ripples
[cohfos1_hpc,cohfos2_hpc]=cellfun(@(equis1,equis2,equis3,equis4,equis5,equis6) co_hfo_spindle(equis1,equis2,equis3,equis4,equis5,equis6),Sx_nhpc,Mx_nhpc,Ex_nhpc,Sx_cortex,Mx_cortex,Ex_cortex,'UniformOutput',false);


Cohfos1_PFC_hpc{k}= ([cohfos1_hpc{:}]);
Cohfos2_PFC_hpc{k}= ([cohfos2_hpc{:}]);
%% Coocur PFC spindle and shuffled HPC ripples.

%Shuffle ripples timestamps
 parfor r=1:1000
[cohfos1_rand_hpc,cohfos2_rand_hpc]=cellfun(@(equis1,equis2,equis3,equis4,equis5,equis6) co_hfo_spindle(equis1,equis2,equis3,equis4,equis5,equis6),Sx_nhpc,Mx_nhpc,Ex_nhpc,Sr_cortex{r},Sr_cortex{r},Er_cortex{r},'UniformOutput',false);

Cohfos1_PFC_hpc_all(k,r)=sum(cellfun('length',cohfos1_rand_hpc));
Cohfos1_PFC_hpc_unique(k,r)=sum(cellfun(@(equis) length(unique(equis)), cohfos1_rand_hpc));
Cohfos2_PFC_hpc_all(k,r)=sum(cellfun('length',cohfos2_rand_hpc));
Cohfos2_PFC_hpc_unique(k,r)=sum(cellfun(@(equis) length(unique(equis)), cohfos2_rand_hpc));

[out_rand_pfc]=coccur_multiplets(cohfos2_rand_hpc);
Out_rand_PFC(k,r,:)=(out_rand_pfc(:,2));
end


%% COOCUR PFC SPINDLE-HPC MULTIPLETS
[out_pfc]=coccur_multiplets(cohfos2_hpc);
Out_PFC{k}=out_pfc;

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
% Find Parietal spindle
[ripple,RipFreq,rip_duration,Mx_par,timeasleep,sig_par,Ex_par,Sx_par,...
  ripple_multiplets_par,RipFreq_multiplets_par,rip_duration_multiplets_par,sig_multiplets_par,Mx_multiplets_par...    
  ]=gui_findspindlesYASA(par,states,yy,multiplets,fn);


si=sig_par(~cellfun('isempty',sig_par));
si=[si{:}];

[x,y,z,~,~,~,l,p]=hfo_specs_spindles(si,timeasleep,fn,0);
%PRE POST ANALYSIS
% [Sx_pre,Mx_pre,Ex_pre,Sx_post,Mx_post,Ex_post]=cellfun(@(equis1,equis2,equis3) pre_post_spindle(equis1,equis2,equis3) ,Sx_hpc,Mx_hpc,Ex_hpc ,'UniformOutput',false);

%%

%Convert signal to 1 sec epochs.
        e_t=1;
        e_samples=e_t*(fn); %fs=1kHz
        ch=length(par);
        nc=floor(ch/e_samples); %Number of epochsw
        NC=[];
        for kk=1:nc
          NC(:,kk)= par(1+e_samples*(kk-1):e_samples*kk);
        end
        vec_bin=states;
        %Use NREM epochs
        vec_bin(vec_bin~=3)=0;
        vec_bin(vec_bin==3)=1;
        %Cluster one values:
        v2=ConsecutiveOnes(vec_bin);
        v_index=find(v2~=0);
        v_values=v2(v2~=0);
        clear v
    for epoch_count=1:length(v_index)
    v{epoch_count,1}=reshape(NC(:, v_index(epoch_count):v_index(epoch_count)+(v_values(1,epoch_count)-1)), [], 1);
    end    
    ti=cellfun(@(equis) reshape(linspace(0, length(equis)-1,length(equis))*(1/fn),[],1) ,v,'UniformOutput',false);

%%
%Find duration of spindles
clear dur_hpc dur_par
for dd=1:length(Mx_par)
        dur_par{dd}=Ex_par{dd}-Sx_par{dd};
    if isempty(Sx_par{dd})
        Sx_par{dd}=nan;
    end

end
dur_par=dur_par.';
% Shuffle timestamps of PAR spindles.
    parfor r=1:1000
        ti_rand=cellfun(@(equis) equis(randperm(size(equis, 1))),ti,'UniformOutput',false);
        Sr_par{r}=cellfun(@(equis1,equis2,equis3) equis3(findclosest(equis1,equis2)).', ti,Sx_par,ti_rand,'UniformOutput',false );
        Er_par{r}=cellfun(@(equis1,equis2,equis3,equis4) equis3(findclosest(equis1,equis2)).'+equis4, ti,Sx_par,ti_rand,dur_par,'UniformOutput',false );
    end

%Change Nans to []
for dd=1:length(Mx_par)
    if isnan(Sx_par{dd})
        Sx_par{dd}=[];
    end
end

%% Find coocurrent PAR spindles and hfos

[cohfos1_g1,cohfos2_g1]=cellfun(@(equis1,equis2,equis3,equis4,equis5,equis6) co_hfo_spindle(equis1,equis2,equis3,equis4,equis5,equis6),Sx_cortex_g1,Mx_cortex_g1,Ex_cortex_g1,Sx_par,Mx_par,Ex_par,'UniformOutput',false);
[cohfos1_g2,cohfos2_g2]=cellfun(@(equis1,equis2,equis3,equis4,equis5,equis6) co_hfo_spindle(equis1,equis2,equis3,equis4,equis5,equis6),Sx_cortex_g2,Mx_cortex_g2,Ex_cortex_g2,Sx_par,Mx_par,Ex_par,'UniformOutput',false);


Cohfos1_PAR_g1{k}=([cohfos1_g1{:}]);
Cohfos2_PAR_g1{k}=([cohfos2_g1{:}]);
Cohfos1_PAR_g2{k}=([cohfos1_g2{:}]);
Cohfos2_PAR_g2{k}=([cohfos2_g2{:}]);

%% Coocur PAR spindle and shuffled hfos.

parfor r=1:1000
[cohfos1_rand_g1,cohfos2_rand_g1]=cellfun(@(equis1,equis2,equis3,equis4,equis5,equis6) co_hfo_spindle(equis1,equis2,equis3,equis4,equis5,equis6),Sx_cortex_g1,Mx_cortex_g1,Ex_cortex_g1,Sr_par{r},Sr_par{r},Er_par{r},'UniformOutput',false);
[cohfos1_rand_g2,cohfos2_rand_g2]=cellfun(@(equis1,equis2,equis3,equis4,equis5,equis6) co_hfo_spindle(equis1,equis2,equis3,equis4,equis5,equis6),Sx_cortex_g2,Mx_cortex_g2,Ex_cortex_g2,Sr_par{r},Sr_par{r},Er_par{r},'UniformOutput',false);


Cohfos1_PAR_all_g1(k,r)=sum(cellfun('length',cohfos1_rand_g1));
Cohfos1_PAR_unique_g1(k,r)=sum(cellfun(@(equis) length(unique(equis)), cohfos1_rand_g1));
Cohfos2_PAR_all_g1(k,r)=sum(cellfun('length',cohfos2_rand_g1));
Cohfos2_PAR_unique_g1(k,r)=sum(cellfun(@(equis) length(unique(equis)), cohfos2_rand_g1));

Cohfos1_PAR_all_g2(k,r)=sum(cellfun('length',cohfos1_rand_g2));
Cohfos1_PAR_unique_g2(k,r)=sum(cellfun(@(equis) length(unique(equis)), cohfos1_rand_g2));
Cohfos2_PAR_all_g2(k,r)=sum(cellfun('length',cohfos2_rand_g2));
Cohfos2_PAR_unique_g2(k,r)=sum(cellfun(@(equis) length(unique(equis)), cohfos2_rand_g2));
end

%% Coocur PAR spindle and HPC ripples
[cohfos1_hpc,cohfos2_hpc]=cellfun(@(equis1,equis2,equis3,equis4,equis5,equis6) co_hfo_spindle(equis1,equis2,equis3,equis4,equis5,equis6),Sx_nhpc,Mx_nhpc,Ex_nhpc,Sx_par,Mx_par,Ex_par,'UniformOutput',false);

Cohfos1_PAR_hpc{k}= ([cohfos1_hpc{:}]);
Cohfos2_PAR_hpc{k}= ([cohfos2_hpc{:}]);

%% Coocur PAR spindle and shuffled HPC ripples.

 parfor r=1:1000
[cohfos1_rand_hpc,cohfos2_rand_hpc]=cellfun(@(equis1,equis2,equis3,equis4,equis5,equis6) co_hfo_spindle(equis1,equis2,equis3,equis4,equis5,equis6),Sx_nhpc,Mx_nhpc,Ex_nhpc,Sr_par{r},Sr_par{r},Er_par{r},'UniformOutput',false);

Cohfos1_PAR_hpc_all(k,r)=sum(cellfun('length',cohfos1_rand_hpc));
Cohfos1_PAR_hpc_unique(k,r)=sum(cellfun(@(equis) length(unique(equis)), cohfos1_rand_hpc));
Cohfos2_PAR_hpc_all(k,r)=sum(cellfun('length',cohfos2_rand_hpc));
Cohfos2_PAR_hpc_unique(k,r)=sum(cellfun(@(equis) length(unique(equis)), cohfos2_rand_hpc));

[out_rand_par]=coccur_multiplets(cohfos2_rand_hpc);
Out_rand_PAR(k,r,:)=(out_rand_par(:,2));
end

%% COOCUR SPINDLE-HPC MULTIPLETS
[out_par]=coccur_multiplets(cohfos2_hpc);
Out_PAR{k}=out_par;


%% Co-occurrent spindles
[cohfos1,cohfos2]=cellfun(@(equis1,equis2,equis3,equis4,equis5,equis6) co_hfo_spindle(equis1,equis2,equis3,equis4,equis5,equis6),Sx_par,Mx_par,Ex_par,Sx_cortex,Mx_cortex,Ex_cortex,'UniformOutput',false);
%Remove repeated values
[cohfos1,cohfos2]=cellfun(@(equis1,equis2) coocur_repeat(equis1,equis2), cohfos1,cohfos2,'UniformOutput',false);

%cohfos1: PAR.
%cohfos2: PFC.

%%
%Co-occurrent spindles
cohf_mx_par=Mx_par(~cellfun('isempty',cohfos1));%Peak values cells where co-occurrent events were found.
cohf_sx_par=Sx_par(~cellfun('isempty',cohfos1));%Peak values cells where co-occurrent events were found.
cohf_ex_par=Ex_par(~cellfun('isempty',cohfos1));%Peak values cells where co-occurrent events were found.

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


%Single Par spindles
v2=cellfun(@(equis1,equis2) single_hfo_get_sample(equis1,equis2),Mx_par,cohfos1,'UniformOutput',false);

Sig_par_single=cellfun(@(equis1,equis2) equis1(equis2),sig_par,v2,'UniformOutput',false);
Sig_par_single=[Sig_par_single{:}];


[single_mx_par_val,single_sx_par_val]=cellfun(@(equis1,equis2,equis3) single_hfos_mx(equis1,equis2,equis3),cohfos1,Mx_par,Sx_par,'UniformOutput',false);
single_mx_par_val=[single_mx_par_val{:}];
single_sx_par_val=[single_sx_par_val{:}];

%%%%
%Cortical COHFOS
cohf_mx_cortex=Mx_cortex(~cellfun('isempty',cohfos2));%Peak values cells where cortex co-occurrent events were found.
cohf_sx_cortex=Sx_cortex(~cellfun('isempty',cohfos2));%Peak values cells where cortex co-occurrent events were found.
cohf_ex_cortex=Ex_cortex(~cellfun('isempty',cohfos2));%Peak values cells where cortex co-occurrent events were found.

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
save('cohfos_spindles_SWR_24.mat','Cohfos2_PAR_hpc_all','Cohfos2_PAR_hpc_unique','Cohfos2_PFC_hpc_all','Cohfos2_PFC_hpc_unique',...
    'Cohfos1_PAR_hpc_all','Cohfos1_PAR_hpc_unique','Cohfos1_PFC_hpc_all','Cohfos1_PFC_hpc_unique',...
       'Cohfos2_PAR_hpc','Cohfos2_PFC_hpc','Cohfos1_PAR_hpc','Cohfos1_PFC_hpc','Out_rand_PAR','Out_rand_PFC','Out_PAR','Out_PFC')

%% Spindle-HPC ripples. random
 TT=table;
    TT.Variables= [mean(Cohfos2_PAR_hpc_all,2).'; mean(Cohfos2_PAR_hpc_unique,2).';...
       mean(Cohfos2_PFC_hpc_all,2).'; mean(Cohfos2_PFC_hpc_unique,2).'; 
    ] ;
    TT.Properties.VariableNames=[cellfun(@(equis) strrep(equis,'_','-'),labelconditions2.','UniformOutput',false).'].';
                writetable(TT,strcat('hpc_ripple_spindles_random','.xls'),'Sheet',1,'Range','A2:L50')    
%% HPC ripples-spindles-multiplets random
avernum_par=mean(Out_rand_PAR,2);
avernum_par=squeeze(avernum_par);
avernum_par=(avernum_par).';
avernum_pfc=mean(Out_rand_PFC,2);
avernum_pfc=squeeze(avernum_pfc);
avernum_pfc=(avernum_pfc).';


TT=table;
TT.Variables=([avernum_par; [NaN NaN NaN NaN] ;avernum_pfc] );
    TT.Properties.VariableNames=[cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false).'].';
                writetable(TT,strcat('hpc_ripple_spindles_multiplets_random','.xls'),'Sheet',1,'Range','A2:L50')    
%% Histograms

% PAR ALL
for n=1:4
  histogram(Cohfos2_PAR_hpc_all(n,:),'FaceColor',[0 0 0])
  hold on
    ylabel('Frequency')
    xlabel('Count')
    Y1 = prctile(Cohfos2_PAR_hpc_all(n,:),5)
    Y2 = prctile(Cohfos2_PAR_hpc_all(n,:),95)

    xline(Y1, '-.k','LineWidth',2)
     xline(Y2, '-.k','LineWidth',2)

  xline(cellfun('length',Cohfos2_PAR_hpc(n)), '-r','LineWidth',2)
pause(1)
       close all

end
%%
% PAR unique
for n=1:4
  histogram(Cohfos2_PAR_hpc_unique(n,:),'FaceColor',[0 0 0])
  hold on
    ylabel('Frequency')
    xlabel('Count')
    Y1 = prctile(Cohfos2_PAR_hpc_unique(n,:),5)
    Y2 = prctile(Cohfos2_PAR_hpc_unique(n,:),95)

    xline(Y1, '-.k','LineWidth',2)
    xline(Y2, '-.k','LineWidth',2)

  xline(cell2mat(cellfun(@(equis) length(unique(equis)),Cohfos2_PAR_hpc(n),'UniformOutput',false)), '-r','LineWidth',2)
     xlim([20 160])
     xticks([20:20:160])
      printing(['Spindle_SWR_coocur_control_PAR_unique_Rat' num2str(Rat) '_' labelconditions2{n}])
      close all

end
%%
% PFC ALL
for n=1:4
  histogram(Cohfos2_PFC_hpc_all(n,:),'FaceColor',[0 0 0])
  hold on
    ylabel('Frequency')
    xlabel('Count')
    Y1 = prctile(Cohfos2_PFC_hpc_all(n,:),5)
    Y2 = prctile(Cohfos2_PFC_hpc_all(n,:),95)

    xline(Y1, '-.k','LineWidth',2)
    xline(Y2, '-.k','LineWidth',2)

  xline(cellfun('length',Cohfos2_PFC_hpc(n)), '-r','LineWidth',2)
     xlim([25 350])
     xticks([25:25: 350])
      printing(['Spindle_SWR_coocur_control_PFC_ALL_Rat' num2str(Rat) '_' labelconditions2{n}])
      close all

end

%%
% PFC unique
for n=1:4
  histogram(Cohfos2_PFC_hpc_unique(n,:),'FaceColor',[0 0 0])
  hold on
    ylabel('Frequency')
    xlabel('Count')
    Y1 = prctile(Cohfos2_PFC_hpc_unique(n,:),5)
    Y2 = prctile(Cohfos2_PFC_hpc_unique(n,:),95)

    xline(Y1, '-.k','LineWidth',2)
    xline(Y2, '-.k','LineWidth',2)

  xline(cell2mat(cellfun(@(equis) length(unique(equis)),Cohfos2_PFC_hpc(n),'UniformOutput',false)), '-r','LineWidth',2)
       xlim([25 175])
       xticks([25:25:175])
      printing(['Spindle_SWR_coocur_control_PFC_unique_Rat' num2str(Rat) '_' labelconditions2{n}])
      close all

end   
  %% Spindle-HFOS
labelconditions2=[    {'plusmaze'}
    {'nl'      }
    {'for'     }
    {'novelty' }];
 TT=table;
    TT.Variables= [mean(Cohfos1_PAR_all_g1,2).'; mean(Cohfos1_PAR_all_g2,2).';...
       mean(Cohfos1_PFC_all_g1,2).'; mean(Cohfos1_PFC_all_g2,2).'; 
    ] ;
    TT.Properties.VariableNames=[cellfun(@(equis) strrep(equis,'_','-'),labelconditions2.','UniformOutput',false).'].';
            writetable(TT,strcat('ripple_spindles_random','.xls'),'Sheet',1,'Range','A2:L10')  

%%
% PAR slow
for n=1:4
  histogram(Cohfos1_PAR_all_g1(n,:),'FaceColor',[0 0 0])
  hold on
    ylabel('Frequency')
    xlabel('Count')
    Y1 = prctile(Cohfos1_PAR_all_g1(n,:),5)
    Y2 = prctile(Cohfos1_PAR_all_g1(n,:),95)

    xline(Y1, '-.k','LineWidth',2)
     xline(Y2, '-.k','LineWidth',2)

  xline(cellfun('length',Cohfos1_PAR_g1(n)), '-r','LineWidth',2)
   xlim([0 4])
   printing(['Spindle_hfo_coocur_control_PAR_g1_Rat' num2str(Rat) '_' labelconditions2{n}])
   close all

end

%% PAR Fast
for n=1:4
  histogram(Cohfos1_PAR_all_g2(n,:),'FaceColor',[0 0 0])
   hold on
  
    ylabel('Frequency')
    xlabel('Count')
    Y1 = prctile(Cohfos1_PAR_all_g2(n,:),5)
    Y2 = prctile(Cohfos1_PAR_all_g2(n,:),95)

    xline(Y1, '-.k','LineWidth',2)
    xline(Y2, '-.k','LineWidth',2)

  xline(cellfun('length',Cohfos1_PAR_g2(n)), '-r','LineWidth',2)
  xlim([0 15])
   printing(['Spindle_hfo_coocur_control_PAR_g2_Rat' num2str(Rat) '_' labelconditions2{n}])
   close all
end
 
%% PFC spindles & slow HFOs
for n=1:4

  histogram(Cohfos1_PFC_all_g1(n,:),'FaceColor',[0 0 0])
  hold on
  
    ylabel('Frequency')
    xlabel('Count')
    Y1 = prctile(Cohfos1_PFC_all_g1(n,:),5)
    Y2 = prctile(Cohfos1_PFC_all_g1(n,:),95)

    xline(Y1, '-.k','LineWidth',2)
    xline(Y2, '-.k','LineWidth',2)

  xline(cellfun('length',Cohfos1_PFC_g1(n)), '-r','LineWidth',2)
  xlim([0 4])
  printing(['Spindle_hfo_coocur_control_PFC_g1_Rat' num2str(Rat) '_' labelconditions2{n}])
  close all
end

%% PFC spindles & fast HFOs
  
for n=1:4
  histogram(Cohfos1_PFC_all_g2(n,:),'FaceColor',[0 0 0])
  hold on
    ylabel('Frequency')
    xlabel('Count')
    Y1 = prctile(Cohfos1_PFC_all_g2(n,:),5)
    Y2 = prctile(Cohfos1_PFC_all_g2(n,:),95)

    xline(Y1, '-.k','LineWidth',2)
    xline(Y2, '-.k','LineWidth',2)

  xline(cellfun('length',Cohfos1_PFC_g2(n)), '-r','LineWidth',2)
  xlim([0 13])
  xticks([0:13])
    printing(['Spindle_hfo_coocur_control_PFC_g2_Rat' num2str(Rat) '_' labelconditions2{n}])
  close all

end
%%
save('cohfos_spindles_hfo_24.mat','Cohfos1_PAR_all_g1','Cohfos1_PAR_all_g2','Cohfos1_PFC_all_g1','Cohfos1_PFC_all_g2',...
       'Cohfos1_PAR_g1','Cohfos1_PAR_g2','Cohfos1_PFC_g1','Cohfos1_PFC_g2')

