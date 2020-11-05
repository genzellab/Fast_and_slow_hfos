%GL_ripples_hfos_control.
%Detects ripples and hfos. Shuffles timestamps to generate a
%null-distribution.
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
tr=getfield(T,strcat('Rat',num2str(Rat)));%Thresholds 
%%
xx={'PAR'};
fn=1000;
%% Get folder names
labelconditions=getfolder;
labelconditions=labelconditions.';
g=labelconditions;

multiplets=[{'singlets'} {'doublets'} {'triplets'} {'quatruplets'} {'pentuplets'} {'sextuplets'} {'septuplets'} {'octuplets'} {'nonuplets'}];
iii=1;
%% Select conditions and sessions.
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

%Load sleep stages files
A = dir('*states*.mat');
A={A.name};

if  ~isempty(A)
       cellfun(@load,A);
else
      error('No Scoring found')    
end

%Find HFOs and create 1000 versions with shuffled timestamps.
[ripple,RipFreq,rip_duration,Mx_cortex,timeasleep,sig_cortex,Ex_cortex,Sx_cortex,...
  ripple_multiplets_cortex,RipFreq_multiplets_cortex,rip_duration_multiplets_cortex,sig_multiplets_cortex,~,Mr ...
  ]=gui_findripples_random(CORTEX,states,xx,tr,multiplets,fn);

si=sig_cortex(~cellfun('isempty',sig_cortex));
si=[si{:}];

All_Par.( strrep(g{k},'-','_'))=si;
[x,y,z,~,~,~,l,p,si_mixed,th]=hfo_specs(si,timeasleep,0,Rat,tr);

%% Cortical HFOs
    hfos_cortex(k)=ripple;
    hfos_cortex_rate(k)=RipFreq;
    hfos_cortex_duration(k)=rip_duration;
    clear ripple RipFreq

%% HPC     
HPC=dir(strcat('*','HPC','*.mat'));
HPC=HPC.name;
HPC=load(HPC);
HPC=getfield(HPC,'HPC');
HPC=HPC.*(0.195);
%Find ripples
[ripple,RipFreq,rip_duration,Mx_hpc,timeasleep,sig_hpc,Ex_hpc,Sx_hpc,...
  ripple_multiplets_hpc,RipFreq_multiplets_hpc,rip_duration_multiplets_hpc,sig_multiplets_hpc,Mx_multiplets_hpc...    
  ]=gui_findripples(HPC,states,{'HPC'},tr,multiplets,fn);


si=sig_hpc(~cellfun('isempty',sig_hpc));
si=[si{:}];

All_HPC.( strrep(g{k},'-','_'))=si;

%% HFC HFOs
hfos_hpc(k)=ripple;
hfos_hpc_rate(k)=RipFreq;
hfos_hpc_duration(k)=rip_duration;

%% Coocurent hfos
[cohfos1,cohfos2]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),Mx_hpc,Mx_cortex,'UniformOutput',false);
%cohfos1: HPC.
%cohfos2: Cortex.
%Common values:
cohfos_count(k)=sum(cellfun('length',cohfos1));
cohfos_rate(k)=sum(cellfun('length',cohfos1))/(timeasleep*(60));


% CONTROL 1000 PERMUTATIONS

for r=1:length(fieldnames(Mr))
[dum_cohfos1,~]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),Mx_hpc,Mr.(['Field_' num2str(r)]),'UniformOutput',false);
random_cohfos_count(k,r)=sum(cellfun('length',dum_cohfos1));
end

%% Random slow and fast Mx timestamp
for i=1:length(fieldnames(Mr)) %1000
    Mx_cortex_g1=Mr.(['Field_' num2str(i)]);
    Mx_cortex_g2=Mr.(['Field_' num2str(i)]);

    row=si_mixed.i1;
    cont=0;
    for ll=1:length(Mx_cortex)

        if ~isempty(Mx_cortex{ll})

            for lll=1:length(Mx_cortex{ll})
                cont=cont+1;
 
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
    Mr_g1.(['Field_' num2str(i)])=Mx_cortex_g1;
    Mr_g2.(['Field_' num2str(i)])=Mx_cortex_g2;
    
end


for r=1:length(fieldnames(Mr))
    [dum_cohfos1_g1,~]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),Mx_hpc,Mr_g1.(['Field_' num2str(r)]),'UniformOutput',false);
    random_cohfos_count_g1(k,r)=sum(cellfun('length',dum_cohfos1_g1));

    [dum_cohfos1_g2,~]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),Mx_hpc,Mr_g2.(['Field_' num2str(r)]),'UniformOutput',false);
    random_cohfos_count_g2(k,r)=sum(cellfun('length',dum_cohfos1_g2));
end


%% Slow and fast hfos
Mx_cortex_g1=Mx_cortex;
Mx_cortex_g2=Mx_cortex;

row=si_mixed.i1;
cont=0;
for ll=1:length(Mx_cortex)

    if ~isempty(Mx_cortex{ll})

        for lll=1:length(Mx_cortex{ll})
            cont=cont+1;

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


%% Coocurrent events
% 
[cohfos1_g1,cohfos2_g1]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),Mx_hpc,Mx_cortex_g1,'UniformOutput',false);
[cohfos1_g2,cohfos2_g2]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),Mx_hpc,Mx_cortex_g2,'UniformOutput',false);

cohfos_count_g1(k)=sum(cellfun('length',cohfos1_g1));
cohfos_rate_g1(k)=sum(cellfun('length',cohfos1_g1))/(timeasleep*(60));

cohfos_count_g2(k)=sum(cellfun('length',cohfos1_g2));
cohfos_rate_g2(k)=sum(cellfun('length',cohfos1_g2))/(timeasleep*(60));

%%

%HPC COHFOS
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

cohf_hpc_dura=cohf_ex_hpc_val-cohf_sx_hpc_val;
cohf_hpc_dura=median(cohf_hpc_dura);
Cohf_hpc_dura(k)=cohf_hpc_dura;

Sig_hpc=sig_hpc(~cellfun('isempty',cohfos1));
Sig_hpc=cellfun(@(equis1,equis2) equis1(equis2),Sig_hpc,coh_samp_hpc,'UniformOutput',false);
Sig_hpc=[Sig_hpc{:}];

%Single HFOs HPC
v2=cellfun(@(equis1,equis2) single_hfo_get_sample(equis1,equis2),Mx_hpc,cohfos1,'UniformOutput',false);

Sig_hpc_single=cellfun(@(equis1,equis2) equis1(equis2),sig_hpc,v2,'UniformOutput',false);
Sig_hpc_single=[Sig_hpc_single{:}];


[single_mx_hpc_val,single_sx_hpc_val]=cellfun(@(equis1,equis2,equis3) single_hfos_mx(equis1,equis2,equis3),cohfos1,Mx_hpc,Sx_hpc,'UniformOutput',false);
single_mx_hpc_val=[single_mx_hpc_val{:}];
single_sx_hpc_val=[single_sx_hpc_val{:}];

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

v2=cellfun(@(equis1,equis2) single_hfo_get_sample(equis1,equis2),Mx_cortex,cohfos2,'UniformOutput',false);

Sig_cortex_single=cellfun(@(equis1,equis2) equis1(equis2),sig_cortex,v2,'UniformOutput',false);
Sig_cortex_single=[Sig_cortex_single{:}];

[single_mx_cortex_val,single_sx_cortex_val]=cellfun(@(equis1,equis2,equis3) single_hfos_mx(equis1,equis2,equis3),cohfos2,Mx_cortex,Sx_cortex,'UniformOutput',false);
single_mx_cortex_val=[single_mx_cortex_val{:}];
single_sx_cortex_val=[single_sx_cortex_val{:}];

progress_bar(k,length(g),f)
    cd ..    
    end

% Histograms of null distribution with actual value indicated. 
%%  
for n=1:4

    histogram(random_cohfos_count(n,:),15)
    ylabel('Frequency')
    xlabel('Count')
    Y1 = prctile(random_cohfos_count(n,:),5)
    Y2 = prctile(random_cohfos_count(n,:),95)
    Y3 = prctile(random_cohfos_count(n,:),2.5)
    Y4 = prctile(random_cohfos_count(n,:),97.5)

    xline(Y1, '-.k','LineWidth',2)
    xline(Y2, '-.k','LineWidth',2)
    xline(Y3, '--k','LineWidth',2)
    xline(Y4, '--k','LineWidth',2)

    xline(cohfos_count(n), '-r','LineWidth',2)

    printing(['Control_Hist' g{n}])
    close all
end

%% Slow PAR
for n=1:4

    histogram(random_cohfos_count_g1(n,:))
    ylabel('Frequency')
    xlabel('Count')
    Y1 = prctile(random_cohfos_count_g1(n,:),5)
    Y2 = prctile(random_cohfos_count_g1(n,:),95)
    Y3 = prctile(random_cohfos_count_g1(n,:),2.5)
    Y4 = prctile(random_cohfos_count_g1(n,:),97.5)

    xline(Y1, '-.k','LineWidth',2)
    xline(Y2, '-.k','LineWidth',2)
    xline(Y3, '--k','LineWidth',2)
    xline(Y4, '--k','LineWidth',2)

    xline(cohfos_count_g1(n), '-r','LineWidth',2)

    printing(['Control_Hist_SLOW_' labelconditions2{n}])
    close all
end

%% Fast PAR
for n=1:4

    histogram(random_cohfos_count_g2(n,:))
    ylabel('Frequency')
    xlabel('Count')
    Y1 = prctile(random_cohfos_count_g2(n,:),5)
    Y2 = prctile(random_cohfos_count_g2(n,:),95)
    Y3 = prctile(random_cohfos_count_g2(n,:),2.5)
    Y4 = prctile(random_cohfos_count_g2(n,:),97.5)

    xline(Y1, '-.k','LineWidth',2)
    xline(Y2, '-.k','LineWidth',2)
    xline(Y3, '--k','LineWidth',2)
    xline(Y4, '--k','LineWidth',2)

    xline(cohfos_count_g2(n), '-r','LineWidth',2)

    printing(['Control_Hist_FAST_' labelconditions2{n}])
    close all
end

%Save values
save('control_random.mat','random_cohfos_count','random_cohfos_count_g1','random_cohfos_count_g2','cohfos_count','cohfos_count_g1','cohfos_count_g2');
%%
vec=random_cohfos_count_g1(:);
 histogram(vec,'FaceColor',[0 0 0])
    ylabel('Frequency')
    xlabel('Count')
    Y1 = prctile(vec,5)
    Y2 = prctile(vec,95)

    xline(Y1, '-.k','LineWidth',2)
     xline(Y2, '-.k','LineWidth',2)

    xline(mean(cohfos_count_g1), '-r','LineWidth',2)
    
y1=prctile([cohfos_count_g1],5)
y2=prctile([cohfos_count_g1],95)

xline(y1, '--r','LineWidth',2)
xline(y2, '--r','LineWidth',2)
%% NORMALIZATION (Z-SCORE-LIKE)

vec_g1_24=random_cohfos_count_g1;
vec_g2_24=random_cohfos_count_g2;
vec_g1_count_24=cohfos_count_g1;
vec_g2_count_24=cohfos_count_g2;

aver_g1_24=vec_g1_24-vec_g1_count_24.';
aver_g1_24=aver_g1_24./std(vec_g1_24.').';

aver_g2_24=vec_g2_24-vec_g2_count_24.';
aver_g2_24=aver_g2_24./std(vec_g2_24.').';

save('aver_norm_24.mat','aver_g1_24','aver_g2_24')
%%
vec_g1_27=random_cohfos_count_g1;
vec_g2_27=random_cohfos_count_g2;
vec_g1_count_27=cohfos_count_g1;
vec_g2_count_27=cohfos_count_g2;

aver_g1_27=vec_g1_27-vec_g1_count_27.';
aver_g1_27=aver_g1_27./std(vec_g1_27.').';

aver_g2_27=vec_g2_27-vec_g2_count_27.';
aver_g2_27=aver_g2_27./std(vec_g2_27.').';

save('aver_norm_27.mat','aver_g1_27','aver_g2_27')
%%
vec_g1_26=random_cohfos_count_g1;
vec_g2_26=random_cohfos_count_g2;
vec_g1_count_26=cohfos_count_g1;
vec_g2_count_26=cohfos_count_g2;

aver_g1_26=vec_g1_26-vec_g1_count_26.';
aver_g1_26=aver_g1_26./std(vec_g1_26.').';

aver_g2_26=vec_g2_26-vec_g2_count_26.';
aver_g2_26=aver_g2_26./std(vec_g2_26.').';

save('aver_norm_26.mat','aver_g1_26','aver_g2_26')

%% Tables with p-values per condition
for n=1:4
vec=aver_g1_24(n,:);
vec2=aver_g2_24(n,:);
pv(n)=(1+sum(vec >=0))/(length(vec)+1);
pv2(n)=(1+sum(vec2 >=0))/(length(vec2)+1);

end
%% p-values for all conditions
vec=aver_g1_24(:,:);
vec=vec(:);
(1+sum(vec >=0))/(length(vec)+1)

vec=aver_g2_24(:,:);
vec=vec(:);
(1+sum(vec >=0))/(length(vec)+1)
%% P-values using all animals

for n=1:4

a_26=aver_g2_26(n,:);
a_27=aver_g2_27(n,:);
a_24=aver_g2_24(n,:);

% a_26=All_norm_post_26(n,:);
% a_27=All_norm_post_27(n,:);
% a_24=All_norm_post_24(n,:);

% a_26=All_norm_post_26(n,:);
% a_27=All_norm_post_27(n,:);
% a_24=All_norm_post_24(n,:);

vec=[a_26 a_27 a_24];

pv(n)=(1+sum(vec >=0))/(length(vec)+1);
end
%% P-values. All animals, all conditions.

a_26=aver_g2_26(:,:);
a_26=a_26(:);
a_27=aver_g2_27(:,:);
a_27=a_27(:);
a_24=aver_g2_24(:,:);
a_24=a_24(:);

vec=[a_26; a_27; a_24]; 
vec=vec(:);
(1+sum(vec >=0))/(length(vec)+1)
