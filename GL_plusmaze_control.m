%GL_plusmaze_control
%Coupling after Plusmaze downsampled to other conditions.
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
xx={'PAR'};
fn=1000;
%% Get folder names
labelconditions=getfolder;
labelconditions=labelconditions.';
g=labelconditions;

multiplets=[{'singlets'} {'doublets'} {'triplets'} {'quatruplets'} {'pentuplets'} {'sextuplets'} {'septuplets'} {'octuplets'} {'nonuplets'}];
iii=1;
%%  Select conditions and sessions.
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

A = dir('*states*.mat');
A={A.name};

if  ~isempty(A)
       cellfun(@load,A);
else
      error('No Scoring found')    
end
%% 
%Find HFOs and create 1000 versions with shuffled timestamps.
[ripple,RipFreq,rip_duration,Mx_cortex,timeasleep,sig_cortex,Ex_cortex,Sx_cortex,...
  ripple_multiplets_cortex,RipFreq_multiplets_cortex,rip_duration_multiplets_cortex,sig_multiplets_cortex,~,Mr ...
  ]=gui_findripples_random(CORTEX,states,xx,tr,multiplets,fn);

si=sig_cortex(~cellfun('isempty',sig_cortex));
si=[si{:}];

All_Par.( strrep(g{k},'-','_'))=si;

%Group actual detections in slow and fast.
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

%% Coocurent events.
[cohfos1,cohfos2]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),Mx_hpc,Mx_cortex,'UniformOutput',false);
%cohfos1: HPC.
%cohfos2: Cortex.
%Common values:
cohfos_count(k)=sum(cellfun('length',cohfos1));
cohfos_rate(k)=sum(cellfun('length',cohfos1))/(timeasleep*(60));

M_cortex.(labelconditions2{k})=Mx_cortex;
M_hpc.(labelconditions2{k})=Mx_hpc;

%% Timestamps of slow and fast HFOs
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

M_cortex_g1.(labelconditions2{k})=Mx_cortex_g1;
M_cortex_g2.(labelconditions2{k})=Mx_cortex_g2;

%% Slow and fast coocurrence with ripple
[cohfos1_g1,cohfos2_g1]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),Mx_hpc,Mx_cortex_g1,'UniformOutput',false);
[cohfos1_g2,cohfos2_g2]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),Mx_hpc,Mx_cortex_g2,'UniformOutput',false);

cohfos_count_g1(k)=sum(cellfun('length',cohfos1_g1));
cohfos_rate_g1(k)=sum(cellfun('length',cohfos1_g1))/(timeasleep*(60));

cohfos_count_g2(k)=sum(cellfun('length',cohfos1_g2));
cohfos_rate_g2(k)=sum(cellfun('length',cohfos1_g2))/(timeasleep*(60));

v2_g1=cellfun(@(equis1,equis2) single_hfo_get_sample(equis1,equis2),Mx_cortex_g1,cohfos2_g1,'UniformOutput',false);
singles_count_g1(k)=sum(cellfun('length',v2_g1));
singles_rate_g1(k)=sum(cellfun('length',v2_g1))/(timeasleep*(60));

v2_g2=cellfun(@(equis1,equis2) single_hfo_get_sample(equis1,equis2),Mx_cortex_g2,cohfos2_g2,'UniformOutput',false);
singles_count_g2(k)=sum(cellfun('length',v2_g2));
singles_rate_g2(k)=sum(cellfun('length',v2_g2))/(timeasleep*(60));


%%


%HPC cooccurrent events
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
  
%% All PAR hfos. 

M=sum(cellfun('length', M_cortex.('plusmaze'))); %Plusmaze

for n=2:4 %Plusmaze vs other conditions.

%Amount of events of condition B.    
N=sum(cellfun('length', M_cortex.(labelconditions2{n}))); %Condition B

for t=1:1000
r=randperm(M,N); %Permute M and take N integers from shuffled M.
cont=0;
M_cortex_dummy=M_cortex;
    for i=1:size(M_cortex.('plusmaze'),1) %EPOCHS
    vec=M_cortex.('plusmaze'){i};

        for j=1:length(vec) %Detections
            cont=cont+1;
            if ~ismember(cont,r)
                M_cortex_dummy.(labelconditions2{1}){i}(j)=NaN;
            end

        end
        
        vec2=M_cortex_dummy.('plusmaze'){i};
     M_cortex_dummy.('plusmaze'){i}=vec2(~isnan(vec2));
    end

    
%Compute cooccurence    
[c1,~]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),M_hpc.('plusmaze'),M_cortex_dummy.('plusmaze'),'UniformOutput',false);

random_count(t)=sum(cellfun('length',c1));
t
end

%Null distribution
histogram(random_count)
    ylabel('Frequency')
    xlabel('Count')
    Y1 = prctile(random_count,5)
    Y2 = prctile(random_count,95)
    Y3 = prctile(random_count,2.5)
    Y4 = prctile(random_count,97.5)
 
    xline(Y1, '-.k','LineWidth',2)
    xline(Y2, '-.k','LineWidth',2)
    xline(Y3, '--k','LineWidth',2)
    xline(Y4, '--k','LineWidth',2)

%Actual value
xline(cohfos_count(n), '-r','LineWidth',2)
 
    printing(['Control_Fed_' labelconditions2{n}])
    close all
end

%% Slow HFOs only.
M=sum(cellfun('length', M_cortex_g1.('plusmaze'))); %Plusmaze

for n=2:4

N=sum(cellfun('length', M_cortex_g1.(labelconditions2{n}))); %Condition B

for t=1:1000
r=randperm(M,N);
cont=0;
M_cortex_g1_dummy=M_cortex_g1;
    for i=1:size(M_cortex_g1.('plusmaze'),1) %EPOCHS
    vec=M_cortex_g1.('plusmaze'){i};

        for j=1:length(vec) %Detections
            cont=cont+1;
            if ~ismember(cont,r)
                M_cortex_g1_dummy.(labelconditions2{1}){i}(j)=NaN;
            end

        end
        
        vec2=M_cortex_g1_dummy.('plusmaze'){i};
     M_cortex_g1_dummy.('plusmaze'){i}=vec2(~isnan(vec2));
    end

    
    
[c1,~]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),M_hpc.('plusmaze'),M_cortex_g1_dummy.('plusmaze'),'UniformOutput',false);

random_count(t)=sum(cellfun('length',c1));
t
end

histogram(random_count)
    ylabel('Frequency')
    xlabel('Count')
    Y1 = prctile(random_count,5)
    Y2 = prctile(random_count,95)
    Y3 = prctile(random_count,2.5)
    Y4 = prctile(random_count,97.5)
% 
    xline(Y1, '-.k','LineWidth',2)
    xline(Y2, '-.k','LineWidth',2)
    xline(Y3, '--k','LineWidth',2)
    xline(Y4, '--k','LineWidth',2)
% 
xline(cohfos_count_g1(n), '-r','LineWidth',2)
% 
    printing(['Control_Fed_Slow_' labelconditions2{n}])
    close all
end
%%
%% FAST HFOs only.
M=sum(cellfun('length', M_cortex_g2.('plusmaze'))); %Plusmaze

for n=2:4

N=sum(cellfun('length', M_cortex_g2.(labelconditions2{n}))); %Condition B

for t=1:1000
r=randperm(M,N);
cont=0;
M_cortex_g2_dummy=M_cortex_g2;
    for i=1:size(M_cortex_g2.('plusmaze'),1) %EPOCHS
    vec=M_cortex_g2.('plusmaze'){i};

        for j=1:length(vec) %Detections
            cont=cont+1;
            if ~ismember(cont,r)
                M_cortex_g2_dummy.(labelconditions2{1}){i}(j)=NaN;
            end

        end
        
        vec2=M_cortex_g2_dummy.('plusmaze'){i};
     M_cortex_g2_dummy.('plusmaze'){i}=vec2(~isnan(vec2));
    end

    
    
[c1,~]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),M_hpc.('plusmaze'),M_cortex_g2_dummy.('plusmaze'),'UniformOutput',false);

random_count(t)=sum(cellfun('length',c1));
t
end

histogram(random_count)
    ylabel('Frequency')
    xlabel('Count')
    Y1 = prctile(random_count,5)
    Y2 = prctile(random_count,95)
    Y3 = prctile(random_count,2.5)
    Y4 = prctile(random_count,97.5)

    xline(Y1, '-.k','LineWidth',2)
    xline(Y2, '-.k','LineWidth',2)
    xline(Y3, '--k','LineWidth',2)
    xline(Y4, '--k','LineWidth',2)
 
xline(cohfos_count_g2(n), '-r','LineWidth',2)
 
    printing(['Control_Fed_Fast_' labelconditions2{n}])
    close all
end
%% Normalized distribution.
%% Slow ripples only.
M=sum(cellfun('length', M_cortex_g1.('plusmaze'))); %Plusmaze

clear aj
for n=2:4
N=sum(cellfun('length', M_cortex_g1.(labelconditions2{n}))); %Condition B
clear random_count
for t=1:1000
r=randperm(M,N);
cont=0;
M_cortex_g1_dummy=M_cortex_g1;
    for i=1:size(M_cortex_g1.('plusmaze'),1) %EPOCHS
    vec=M_cortex_g1.('plusmaze'){i};

        for j=1:length(vec) %Detections
            cont=cont+1;
            if ~ismember(cont,r)
                M_cortex_g1_dummy.(labelconditions2{1}){i}(j)=NaN;
            end

        end
        
        vec2=M_cortex_g1_dummy.('plusmaze'){i};
     M_cortex_g1_dummy.('plusmaze'){i}=vec2(~isnan(vec2));
    end

    
%Compute cooccurrence    
[c1,~]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),M_hpc.('plusmaze'),M_cortex_g1_dummy.('plusmaze'),'UniformOutput',false);

random_count(t)=sum(cellfun('length',c1));
t
end

clear r_count
%Normalization
r_count=(random_count-cohfos_count_g1(n))/std(random_count);
Fed_slow_24(n-1,:)=r_count;
set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultTextFontName', 'Arial')

histogram(r_count,[-10.5:1:10.5],'FaceColor',[0 0 0])
    ylabel('Frequency','FontSize',10,'FontName','Arial')
    xlabel('Count','FontSize',10,'FontName','Arial')
    Y1 = prctile(r_count,5)
    Y2 = prctile(r_count,95)

    xline(Y1, '-.k','LineWidth',2)
    xline(Y2, '-.k','LineWidth',2)

xline(0, '-r','LineWidth',2)
xlim([-10 10])
xticks([-10:2:10])

aj(n)=(1+sum(r_count >=0))/(length(r_count)+1)

set(gca,'FontName','Arial')
set(gca,'FontSize',10)
ax = gca;
ax.YAxis.FontSize = 18 %for y-axis 
ax.XAxis.FontSize = 18 %for y-axis 
ax.XAxis.FontName='Arial';
ax.XAxis.FontName='Arial';

     printing(['Control_Fed_Slow_' num2str(Rat) '_' labelconditions2{n}])
    close all
end

%%
%% FAST ripples only.
M=sum(cellfun('length', M_cortex_g2.('plusmaze'))); %Plusmaze

clear aj
for n=2:4
N=sum(cellfun('length', M_cortex_g2.(labelconditions2{n}))); %Condition B
clear random_count
for t=1:1000
r=randperm(M,N);
cont=0;
M_cortex_g2_dummy=M_cortex_g2;
    for i=1:size(M_cortex_g2.('plusmaze'),1) %EPOCHS
    vec=M_cortex_g2.('plusmaze'){i};

        for j=1:length(vec) %Detections
            cont=cont+1;
            if ~ismember(cont,r)
                M_cortex_g2_dummy.(labelconditions2{1}){i}(j)=NaN;
            end

        end
        
        vec2=M_cortex_g2_dummy.('plusmaze'){i};
     M_cortex_g2_dummy.('plusmaze'){i}=vec2(~isnan(vec2));
    end

    
    
[c1,~]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),M_hpc.('plusmaze'),M_cortex_g2_dummy.('plusmaze'),'UniformOutput',false);

random_count(t)=sum(cellfun('length',c1));
t
end
clear r_count
r_count=(random_count-cohfos_count_g2(n))/std(random_count);
Fed_Fast_24(n-1,:)=r_count;
set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultTextFontName', 'Arial')

histogram(r_count,[-10.5:1:10.5],'FaceColor',[0 0 0])
    ylabel('Frequency','FontSize',10,'FontName','Arial')
    xlabel('Count','FontSize',10,'FontName','Arial')
    Y1 = prctile(r_count,5)
    Y2 = prctile(r_count,95)
    xline(Y1, '-.k','LineWidth',2)
    xline(Y2, '-.k','LineWidth',2)
xline(0, '-r','LineWidth',2)
xlim([-10 10])
xticks([-10:2:10])

aj(n)=(1+sum(r_count >=0))/(length(r_count)+1)

set(gca,'FontName','Arial')
set(gca,'FontSize',10)
ax = gca;
ax.YAxis.FontSize = 18 %for y-axis 
ax.XAxis.FontSize = 18 %for y-axis 
ax.XAxis.FontName='Arial';
ax.XAxis.FontName='Arial';

     printing(['Control_Fed_Fast_' num2str(Rat) '_' labelconditions2{n}])
    close all
end

%% Compute p-values
vec=Fed_Fast_24;
vec=vec(:);

(1+sum(vec >=0))/(length(vec)+1)
%%

for n=1:3
a_26=Fed_Fast_26(n,:);
a_27=Fed_Fast_27(n,:);
a_24=Fed_Fast_24(n,:);

vec=[a_26 a_27 a_24];

aj(n)=(1+sum(vec >=0))/(length(vec)+1);
end
%%
a_26=Fed_Fast_26(:,:);
a_26=a_26(:);
a_27=Fed_Fast_27(:,:);
a_27=a_27(:);
a_24=Fed_Fast_24(:,:);
a_24=a_24(:);

vec=[a_26; a_27; a_24]; %a_26; a_27; a_24
vec=vec(:);
(1+sum(vec >=0))/(length(vec)+1)
