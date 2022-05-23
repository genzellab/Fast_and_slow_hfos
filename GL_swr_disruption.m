%gui_threshold_ripples
% cd('/home/adrian/Dropbox/new_plusmaze')
% load('plusmaze_21_25.mat')
cd('/home/adrian/Documents/GitHub/CorticoHippocampal/Fast_and_slow_hfos')
load('load_me_first.mat')
%% Find location
close all
% dname=uigetdir([],'Select folder with Matlab data containing all rats.');
% cd(dname)
%cd('/home/adrian/Dropbox/jukebox/Desktop/SWRD_extra_rats/downsampled_files')
cd('/media/adrian/6aa1794c-0320-4096-a7df-00ab0ba946dc/jukebox/Desktop/SWRD_extra_rats/downsampled_files')
%%
% z= zeros(length(label1),length(rats));
% [T]=gui_table_channels(z,rats,label1,'Threholds');

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
        % Ask for brain area.
% xx = inputdlg({'Cortical Brain area'},...
%               'Type your selection', [1 30]); 
xx={'PAR'};
fn=1000;
%%
gg=getfolder;
gg=gg.';
if size(label1,1)~=3  % IF not Plusmaze
    gg(ismember(gg,'OR_N'))=[];
    gg(ismember(gg,'OD_N(incomplete)'))=[];
    gg=sort(gg); %Sort alphabetically.
    labelconditions2=gg;
    gg(ismember(gg,'CN'))={'CON'};
end
labelconditions=gg;

%% Select experiment to perform. 
inter=1;
%Select length of window in seconds:
ro=[1200];
coher=0;
selectripples=1;
notch=0; %Might need to be 1.
nrem=3;
level=1;

multiplets=[{'singlets'} {'doublets'} {'triplets'} {'quatruplets'} {'pentuplets'} {'sextuplets'} {'septuplets'} {'octuplets'} {'nonuplets'}];
%%
iii=1;
%for iii=1:length(labelconditions) 

    if size(label1,1)~=3  % IF not Plusmaze

        cd( labelconditions2{iii})
        g=getfolder;

        if iii==1
            answer = questdlg('Should we use all trials?', ...
                'Trial selection', ...
                'Use all','Select trials','Select trials');

            % Handle response
            switch answer
                case 'Use all'
                    disp(['Using all.'])
                    an=[];
                case 'Select trials'
                    prompt = {['Enter trials name common word without index:' sprintf('\n') '(Use commas for multiple names)']};
                    dlgtitle = 'Input';
                    dims = [1 35];
                    %definput = {'20','hsv'};
                    an = inputdlg(prompt,dlgtitle,dims);
                    %an=char(an);
            %        g=g(contains(g,{'PT'}));
            end

        end

        if ~isempty(an)
        g=g(contains(g,strsplit(an{1},',')));
        end
  
    else
      g=gg;  
    end
  %% Colormap
       % n=length(g);
myColorMap=jet(length(g));              
%%
f=waitbar(0,'Please wait...');
    for k=1:length(g)
        cd(g{k})
%(level,nrem,notch,w,lepoch)
A = dir('*states*.mat');
A={A.name};

if  ~isempty(A)
       cellfun(@load,A);
else
      error('No Scoring found')    
end
%% HPC     
HPC=dir(strcat('*','HPC','*.mat'));
HPC=HPC.name;
HPC=load(HPC);
%HPC=HPC.HPC;
HPC=getfield(HPC,'HPC');
HPC=HPC.*(0.195);
%xo
 %Remove 50Hz and harmonics
[HPC]=notch_filter(HPC,fn);
%xo 

[ripple,RipFreq,rip_duration,Mx_hpc,timeasleep,sig_hpc,Ex_hpc,Sx_hpc,...
  ripple_multiplets_hpc,RipFreq_multiplets_hpc,rip_duration_multiplets_hpc,sig_multiplets_hpc,Mx_multiplets_hpc...    
  V,Mono,wa,wa2]=gui_findripples_swrd(HPC,states,{'HPC'},tr,multiplets,fn,Rat);



%CORTEX,states,xx,tr,multiplets,fn

si=sig_hpc(~cellfun('isempty',sig_hpc));
si=[si{:}];

% plot_hfo(si,Mx_hpc,Sx_hpc,label1{1})
% title(['HFO HPC  ' strrep(g{k},'_','-')])
% cd ..
% printing(['HFO HPC  ' strrep(g{k},'_','-')])
% close all
% cd(g{k})
% xo
All_HPC.( strrep(g{k},'-','_'))=si;
[x,y,z,~,~,~,l,p]=hfo_specs_hpc(si,timeasleep,0);
% cd ..
% % printing(['SWRD_Histograms_HPC_Probability_' g{k}]);
% close all
% cd(g{k})

fi_hpc(k)=x;
fa_hpc(k)=y;
amp_hpc(k)=z;
auc_hpc(k)=l;
p2p_hpc(k)=p;
% %Instantaneous frequency.
% x=cellfun(@(equis) mean(instfreq(equis,1000)) ,si,'UniformOutput',false);
% x=cell2mat(x);
% x=median(x);
% fi_hpc(k)=x;
% %Average frequency
% y=cellfun(@(equis) (meanfreq(equis,1000)) ,si,'UniformOutput',false);
% y=cell2mat(y);
% y=median(y);
% fa_hpc(k)=y;
% 
% %Amplitude
% z=cellfun(@(equis) max(abs(hilbert(equis))) ,si,'UniformOutput',false);
% z=cell2mat(z);
% z=median(z);
% amp_hpc(k)=z;

% Mx_cortex(~cellfun('isempty',Mx_cortex))
%% HPC ripples
hfos_hpc(k)=ripple;
hfos_hpc_rate(k)=RipFreq;
hfos_hpc_duration(k)=rip_duration;

%Multiplets    
for ll=1:length(multiplets)
   eval(['hfos_hpc_' multiplets{ll} '(k)=ripple_multiplets_hpc.' multiplets{ll} ';']) 
   eval(['hfos_hpc_rate_' multiplets{ll} '(k)=RipFreq_multiplets_hpc.' multiplets{ll} ';']) 
   eval(['hfos_hpc_duration_' multiplets{ll} '(k)=rip_duration_multiplets_hpc.' multiplets{ll} ';'])    
end
%xo

%%
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
CORTEX=CORTEX.*(0.195);


 %Remove 50Hz and harmonics
[CORTEX]=notch_filter(CORTEX,fn);
 xo 
[si,Sx,Ex,Mx,ti,Mono_cortex,V_cortex,timeasleep]=gui_findripples_swrd_part1(CORTEX,states,xx,tr,fn,wa2,Rat);

[~,~,~,~,~,~,~,~,~,th,PCA_features]=hfo_specs(si,timeasleep,0,Rat,tr);
%Remove frequency features to avoid bias.
 %PCA_features=PCA_features(:,3:end);
 PCA_features=PCA_features(:,[3 4 6 7]);

rng('Default')
[coeff,score,latent,tsquared,explainedVar] = pca(PCA_features);
% bar(explainedVar)
% title('Explained Variance: More than 90% explained by first two principal components')
% ylabel('PC')
% %%
% imagesc(coeff)

%xo
% %%
% % Retain first two principal components
% yeastPC = score(:,1:2);
% 
% figure()
% rng('default')
% % rng(3)
% [clusters, centroid] = kmeans(yeastPC,2);
% 
% gscatter(yeastPC(:,1),yeastPC(:,2),clusters)
% legend('location','southeast')
% xlabel('First Principal Component');
% ylabel('Second Principal Component');
% title('Principal Component Scatter Plot with Colored Clusters');
%% DB scan approach
yeastPC = score(:,1:2);
 rng('Default')
%rng(1)

idx = dbscan(yeastPC,10,10);
gscatter(yeastPC(:,1),yeastPC(:,2),idx)
legend('location','southeast')
xlabel('First Principal Component');
ylabel('Second Principal Component');
title('Principal Component Scatter Plot with Colored Clusters');
close all
% %% DB scan approach (Z-scored)
% aver=zscore(yeastPC);
% idx = dbscan(aver,.25,10); %.25, %10
% gscatter(aver(:,1),aver(:,2),idx)
% legend('location','southeast')
% xlabel('First Principal Component');
% ylabel('Second Principal Component');
% title('Principal Component Scatter Plot with Colored Clusters');

%% Find lower left corner
clustercoef=unique(idx);
clustercoef=clustercoef(clustercoef~=-1);

KL=[100 1000];
if ~isempty(clustercoef)
    for kl=1:length(clustercoef)

    if mean(yeastPC(idx==clustercoef(kl),1))<KL(1) %& mean(yeastPC(idx==clustercoef(kl),1))<KL(1)
        cl=clustercoef(kl); %Cluster ID at the lower left corner.
        KL=[mean(yeastPC(idx==clustercoef(kl),1))  mean(yeastPC(idx==clustercoef(kl),2))];
    end

    end
else
   cl=1;
   idx=idx+2;
end

sum(idx==cl)
% %% REMOVE all consecutive detections
% averno=idx;
% averno(averno~=cl)=0;
% averno=averno/cl;
% avern=ConsecutiveOnes(averno);
% 
% % find(avern>floor(median(avern(avern~=0))))
% % idx(find(avern>1))=-1;
% nx=find(avern>floor(median(avern(avern~=0))));
% Vec_ind=[];
% for hj=1:length(nx)
%     Vec_ind=[Vec_ind nx(hj):nx(hj)+avern(nx(hj))-1];
%     
% end
% idx(Vec_ind)=-1;
% sum(idx==cl)
%% Take some of the consecutive detections
%xo
averno=idx;
averno(averno~=cl)=0;
averno=averno/cl;
avern=ConsecutiveOnes(averno);
nx=find(avern>=min(avern(avern~=0)));

Vec_ind=[];
for hj=1:length(nx)
    Vec_ind=[Vec_ind nx(hj)+min(avern(avern~=0)):nx(hj)+avern(nx(hj))-1];
    
end
idx(Vec_ind)=-1;
sum(idx==cl)


% %% 3D
% yeastPC = score(:,1:3);
% 
% figure()
% 
% scatter3(yeastPC(:,1),yeastPC(:,2),yeastPC(:,3))
% legend('location','southeast')
% xlabel('First Principal Component');
% ylabel('Second Principal Component');
% zlabel('Third Principal Component');
% 
% title('Principal Component Scatter Plot with Colored Clusters');
% %%
% [sum((score(:,1)<0))
% sum((score(:,1)<0).*(score(:,2)<0))
% 
% sum((score(:,1)<0).*(score(:,2)<0).*(score(:,3)<0))
% 
% sum((score(:,1)<0).*(score(:,2)<0).*(score(:,3)<0).*(score(:,4)<0))]
% %%
% [sum((score(:,1)< prctile(score(:,1),10)))
% sum((score(:,1)< prctile(score(:,1),20)))
% sum((score(:,1)< prctile(score(:,1),30)))
% sum((score(:,1)< prctile(score(:,1),40)))
% sum((score(:,1)< prctile(score(:,1),50)))
% 
% ]
% %%
% [
%     sum( zscore(score(:,1))<0 )
% sum( (zscore(score(:,1))<0).*(zscore(score(:,2))<0)   )
% 
% sum( (zscore(score(:,1))<0).*(zscore(score(:,2))<0).*(zscore(score(:,3))<0)     )
% 
% sum( (zscore(score(:,1))<0).*(zscore(score(:,2))<0).*(zscore(score(:,3))<0) .*(zscore(score(:,4))<0)     )
% ]
%%

% sum((score(:,1)<0).*(score(:,2)<0))
% 
% sum((score(:,1)<0).*(score(:,2)<0).*(score(:,3)<0))
% 
% sum((score(:,1)<0).*(score(:,2)<0).*(score(:,3)<0).*(score(:,4)<0))]

%%
% find((yeastPC(:,1)<0).*(yeastPC(:,2)<0).*(yeastPC(:,3)<0))

%% condition
%sum(yeastPC(:,1)<=0)/length(yeastPC(:,1)) %60
% 1.5*std(yeastPC(:,1))
%%
% SI=si(find(clusters==1).')
% %%
% skewness(yeastPC(:,2))
SI=si(idx==cl);
%% Correction
% av=(cellfun(@(equis) max(abs(equis)),si));
% index_false_positives=find(isoutlier(av)); %This index is with respect to si, which is the number of events, not epochs.
index_false_positives=find(idx~=cl);
%xo
%%
[ripple,RipFreq,rip_duration,Mx_cortex,timeasleep,sig_cortex,Ex_cortex,Sx_cortex,...
  ripple_multiplets_cortex,RipFreq_multiplets_cortex,rip_duration_multiplets_cortex,sig_multiplets_cortex,~ ...
  ]=gui_findripples_swrd_part2(Sx,Ex,Mx,index_false_positives,ti,Mono_cortex,multiplets,timeasleep);
%% Remove consecutive large stimulation peaks
ssi=si;
for chale=1:length(si)
    [a,b]=findpeaks(movmad(abs(si{chale}),10),'MinPeakDistance',15','MinPeakHeight',1.5);
    if length(a)<3
        ssi(chale)={NaN};
    end
end
%%

SSI=cellfun(@(equis) length(equis)  ,ssi,'UniformOutput',false);
a=[SSI{:}];
A=a~=1;
index_false_positives2=find(A).';

[ripple,RipFreq,rip_duration,Mx_cortex,timeasleep,sig_cortex,Ex_cortex,Sx_cortex,...
  ripple_multiplets_cortex,RipFreq_multiplets_cortex,rip_duration_multiplets_cortex,sig_multiplets_cortex,~ ...
  ]=gui_findripples_swrd_part2(Sx_cortex,Ex_cortex,Mx_cortex,index_false_positives2,ti,Mono_cortex,multiplets,timeasleep);


%%
% [ripple,RipFreq,rip_duration,Mx_cortex,timeasleep,sig_cortex,Ex_cortex,Sx_cortex,...
%   ripple_multiplets_cortex,RipFreq_multiplets_cortex,rip_duration_multiplets_cortex,sig_multiplets_cortex,~, ...
%   V,Mono]=gui_findripples_swrd(CORTEX,states,xx,tr,multiplets,fn);

si=sig_cortex(~cellfun('isempty',sig_cortex));
si=[si{:}];
%xo
% plot_hfo(si,Mx_cortex,Sx_cortex,label1{2})
% title(['HFO Cortex  ' strrep(g{k},'_','-')])
% cd ..
% printing(['HFO Cortex  ' strrep(g{k},'_','-')])
% close all
% cd(g{k})
All_Par.( strrep(g{k},'-','_'))=si;
% All_timeasleep.( strrep(g{k},'-','_'))=timeasleep;
%xo
[x,y,z,~,~,~,l,p,si_mixed,th]=hfo_specs(si,timeasleep,0,Rat,tr);
% cd ..
% % printing(['SWRD_Histograms_Cortex_Count_' g{k}]);
% close all
% cd(g{k})

fi_cortex(k)=x;
fa_cortex(k)=y;
amp_cortex(k)=z;
auc_cortex(k)=l;
p2p_cortex(k)=p;
%xo
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
    
    
%     C = cellfun(@minus,Ex_pfc,Sx_pfc,'UniformOutput',false);
%     CC=([C{:}]);
%     hfos_pfc_duration(k)=median(CC);
%% Coocurent hfos
[cohfos1,cohfos2]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),Mx_hpc,Mx_cortex,'UniformOutput',false);
%cohfos1: HPC.
%cohfos2: Cortex.
%Common values:
cohfos_count(k)=sum(cellfun('length',cohfos1));
cohfos_rate(k)=sum(cellfun('length',cohfos1))/(timeasleep*(60));

%Multiplet cohfos
for ll=1:length(multiplets)
[cohfos1_multiplets.(multiplets{ll}),cohfos2_multiplets.(multiplets{ll})]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),Mx_multiplets_hpc.(multiplets{ll}).',Mx_cortex,'UniformOutput',false);
cohfos_count_multiplets.(multiplets{ll})(k)=sum(cellfun('length',cohfos1_multiplets.(multiplets{ll})));
cohfos_rate_multiplets.(multiplets{ll})(k)=sum(cellfun('length',cohfos1_multiplets.(multiplets{ll})))/(timeasleep*(60));
end

%% Mixed distribution (Average freq) coHFOs
Mx_cortex_g1=Mx_cortex;
Mx_cortex_g2=Mx_cortex;

row=si_mixed.i1;
cont=0;
for ll=1:length(Mx_cortex)
% cont=cont+length(Mx_cortex{ll});

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

[cohfos1_g1,cohfos2_g1]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),Mx_hpc,Mx_cortex_g1,'UniformOutput',false);
[cohfos1_g2,cohfos2_g2]=cellfun(@(equis1,equis2) co_hfo(equis1,equis2),Mx_hpc,Mx_cortex_g2,'UniformOutput',false);

cohfos_count_g1(k)=sum(cellfun('length',cohfos1_g1));
cohfos_rate_g1(k)=sum(cellfun('length',cohfos1_g1))/(timeasleep*(60));

cohfos_count_g2(k)=sum(cellfun('length',cohfos1_g2));
cohfos_rate_g2(k)=sum(cellfun('length',cohfos1_g2))/(timeasleep*(60));
%xo

v2_g1=cellfun(@(equis1,equis2) single_hfo_get_sample(equis1,equis2),Mx_cortex_g1,cohfos2_g1,'UniformOutput',false);
singles_count_g1(k)=sum(cellfun('length',v2_g1));
singles_rate_g1(k)=sum(cellfun('length',v2_g1))/(timeasleep*(60));


v2_g2=cellfun(@(equis1,equis2) single_hfo_get_sample(equis1,equis2),Mx_cortex_g2,cohfos2_g2,'UniformOutput',false);
singles_count_g2(k)=sum(cellfun('length',v2_g2));
singles_rate_g2(k)=sum(cellfun('length',v2_g2))/(timeasleep*(60));

all_count_g1(k)=[sum(cellfun('length',v2_g1))+sum(cellfun('length',cohfos1_g1))];
all_rate_g1(k)=[sum(cellfun('length',v2_g1))+sum(cellfun('length',cohfos1_g1))]/(timeasleep*(60));

all_count_g2(k)=[sum(cellfun('length',v2_g2))+sum(cellfun('length',cohfos1_g2))];
all_rate_g2(k)=[sum(cellfun('length',v2_g2))+sum(cellfun('length',cohfos1_g2))]/(timeasleep*(60));

all_count(k)=[sum(cellfun('length',v2_g1))+sum(cellfun('length',cohfos1_g1))]+[sum(cellfun('length',v2_g2))+sum(cellfun('length',cohfos1_g2))];
all_rate(k)=all_count(k)/(timeasleep*(60));
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
%xo
Sig_hpc=sig_hpc(~cellfun('isempty',cohfos1));
Sig_hpc=cellfun(@(equis1,equis2) equis1(equis2),Sig_hpc,coh_samp_hpc,'UniformOutput',false);
Sig_hpc=[Sig_hpc{:}];
%xo

% plot_hfo(Sig_hpc,{cohf_mx_hpc_val},{cohf_sx_hpc_val},label1{1})
% title(['coHFO HPC envelope  ' strrep(g{k},'_','-')])
% cd ..
% printing(['coHFO HPC envelope ' strrep(g{k},'_','-')])
% close all
% cd(g{k})



[x,y,z,w,h,q,l,p]=hfo_specs_hpc(Sig_hpc,timeasleep,0);
fi_cohfo_hpc(k)=x;
fa_cohfo_hpc(k)=y;
amp_cohfo_hpc(k)=z;
count_cohfo_hpc(k)=w;
rate_cohfo_hpc(k)=h;
dura_cohfo_hpc(k)=q;
auc_cohfo_hpc(k)=l;
p2p_cohfo_hpc(k)=p;
%Single HFOs HPC
%[v2]=single_hfo_get_sample(Mx_hpc{1},cohfos1{1});
v2=cellfun(@(equis1,equis2) single_hfo_get_sample(equis1,equis2),Mx_hpc,cohfos1,'UniformOutput',false);

Sig_hpc_single=cellfun(@(equis1,equis2) equis1(equis2),sig_hpc,v2,'UniformOutput',false);
Sig_hpc_single=[Sig_hpc_single{:}];


[single_mx_hpc_val,single_sx_hpc_val]=cellfun(@(equis1,equis2,equis3) single_hfos_mx(equis1,equis2,equis3),cohfos1,Mx_hpc,Sx_hpc,'UniformOutput',false);
single_mx_hpc_val=[single_mx_hpc_val{:}];
single_sx_hpc_val=[single_sx_hpc_val{:}];
% xo


% plot_hfo(Sig_hpc_single,{single_mx_hpc_val},{single_sx_hpc_val},label1{1})
% title(['Single HPC envelope ' strrep(g{k},'_','-')])
% cd ..
% printing(['Single HPC envelope ' strrep(g{k},'_','-')])
% close all
% cd(g{k})

[x,y,z,w,h,q,l,p]=hfo_specs_hpc(Sig_hpc_single,timeasleep,0);
fi_single_hpc(k)=x;
fa_single_hpc(k)=y;
amp_single_hpc(k)=z;
count_single_hpc(k)=w;
rate_single_hpc(k)=h;
dura_single_hpc(k)=q;
auc_single_hpc(k)=l;
p2p_single_hpc(k)=p;


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
% xo

% plot_hfo(Sig_cortex,{cohf_mx_cortex_val},{cohf_sx_cortex_val},label1{2})
% title(['coHFO cortex envelope ' strrep(g{k},'_','-')])
% cd ..
% printing(['coHFO cortex envelope ' strrep(g{k},'_','-')])
% close all
% cd(g{k})


[x,y,z,w,h,q,l,p]=hfo_specs(Sig_cortex,timeasleep,0,Rat,tr);
fi_cohfo_cortex(k)=x;
fa_cohfo_cortex(k)=y;
amp_cohfo_cortex(k)=z;
count_cohfo_cortex(k)=w;
rate_cohfo_cortex(k)=h;
dura_cohfo_cortex(k)=q;
auc_cohfo_cortex(k)=l;
p2p_cohfo_cortex(k)=p;
%Single HFOs Cortex
%[v2]=single_hfo_get_sample(Mx_hpc{1},cohfos1{1});
v2=cellfun(@(equis1,equis2) single_hfo_get_sample(equis1,equis2),Mx_cortex,cohfos2,'UniformOutput',false);

Sig_cortex_single=cellfun(@(equis1,equis2) equis1(equis2),sig_cortex,v2,'UniformOutput',false);
Sig_cortex_single=[Sig_cortex_single{:}];
%xo
[single_mx_cortex_val,single_sx_cortex_val]=cellfun(@(equis1,equis2,equis3) single_hfos_mx(equis1,equis2,equis3),cohfos2,Mx_cortex,Sx_cortex,'UniformOutput',false);
single_mx_cortex_val=[single_mx_cortex_val{:}];
single_sx_cortex_val=[single_sx_cortex_val{:}];



% plot_hfo(Sig_cortex_single,{single_mx_cortex_val},{single_sx_cortex_val},label1{2})
% title(['Single cortex envelope ' strrep(g{k},'_','-')])
% cd ..
% printing(['Single cortex envelope ' strrep(g{k},'_','-')])
% close all
% cd(g{k})

[x,y,z,w,h,q,l,p]=hfo_specs(Sig_cortex_single,timeasleep,0,Rat,tr);
fi_single_cortex(k)=x;
fa_single_cortex(k)=y;
amp_single_cortex(k)=z;
count_single_cortex(k)=w;
rate_single_cortex(k)=h;
dura_single_cortex(k)=q;
auc_single_cortex(k)=l;
p2p_single_cortex(k)=p;
%%
% xo
progress_bar(k,length(g),f)
    cd ..    
    end
  xo
% All_timeasleep
% All_Par
% A_cell = struct2cell(All_Par);
% All_Par_24_35=[A_cell{:}];
% %%
% All_40=[All_Par_26 All_Par_27 All_Par_24_40];
% All_35=[All_Par_26 All_Par_27 All_Par_24_35];
% %%
% 
% si=All_40;
% hfo_specs(si,1,1)
% %printing('Histogram_All_Par_40_probability')
% printing('Histogram_All_Par_40_count')
% 
% %%
% si=All_35;
% hfo_specs(si,1,1)
% % printing('Histogram_All_Par_35_probability')
% printing('Histogram_All_Par_35_count')

% %AUC
% TT=table;
% TT.Variables=    [[{'All_cortex'};{'All_HPC'};{'Cohfo_cortex'};{'Cohfo_hpc'};{'Single_cortex'};{'Single_HPC'}] num2cell([auc_cortex;auc_hpc;auc_cohfo_cortex;auc_cohfo_hpc;auc_single_cortex;auc_single_hpc])];
% TT.Properties.VariableNames=['HFO Type';cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)].';
% writetable(TT,strcat('AUC_',num2str(tr(1)),'_',num2str(tr(2)),'.xls'),'Sheet',1,'Range','A2:L10')    
% 
% %Peak-to-Peak distance
% TT=table;
% TT.Variables=    [[{'All_cortex'};{'All_HPC'};{'Cohfo_cortex'};{'Cohfo_hpc'};{'Single_cortex'};{'Single_HPC'}] num2cell([p2p_cortex;p2p_hpc;p2p_cohfo_cortex;p2p_cohfo_hpc;p2p_single_cortex;p2p_single_hpc])];
% TT.Properties.VariableNames=['HFO Type';cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)].';
% writetable(TT,strcat('P2P_',num2str(tr(1)),'_',num2str(tr(2)),'.xls'),'Sheet',1,'Range','A2:L10')    
% %xo


%Cortex
% c = categorical(cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)); 
% bar(c,hfos_cortex)
% ylabel('Number of HFOs')
% title(xx{1})
% 
%     if size(label1,1)~=3  % IF not Plusmaze 
%       string=strcat('HFOs_counts_',xx{1},'_Rat',num2str(Rat),'_',labelconditions{iii}); 
%     else
% %         if strcmp(xx{1},'HPC')
% %                   string=strcat('HFOs_counts_',xx{1},'_Rat',num2str(Rat),'_',num2str(tr(1)));         
% %         else
%                   string=strcat('HFOs_counts_',xx{1},'_Rat',num2str(Rat),'_',num2str(tr(2)));         
% %         end
%     end
% 
%     printing(string)
%     close all
%rate
c = categorical(cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)); 
bar(c,hfos_cortex_rate)
ylabel('HFOs per second')
title(xx{1})


    if size(label1,1)~=3  % IF not Plusmaze 
      string=strcat('HFOs_rate_',xx{1},'_Rat',num2str(Rat),'_',labelconditions{iii}); 
    else
%         if strcmp(xx{1},'HPC')
%                   string=strcat('HFOs_rate_',xx{1},'_Rat',num2str(Rat),'_',num2str(tr(1)));         
%         else
                  string=strcat('HFOs_rate_',xx{1},'_Rat',num2str(Rat),'_',num2str(tr(2)));         
%         end
        
    end

    printing(string)
    close all
    
%Average Frequency
% c = categorical(cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)); 
% bar(c,fa_cortex)
% ylabel('Average frequency')
% title(xx{1})
% ylim([100 300])
% 
%     if size(label1,1)~=3  % IF not Plusmaze 
%       string=strcat('HFOs_counts_',xx{1},'_Rat',num2str(Rat),'_',labelconditions{iii}); 
%     else
% %         if strcmp(xx{1},'HPC')
% %                   string=strcat('HFOs_counts_',xx{1},'_Rat',num2str(Rat),'_',num2str(tr(1)));         
% %         else
%                   string=strcat('HFOs_average_frequency_',xx{1},'_Rat',num2str(Rat),'_',num2str(tr(2)));         
% %         end
%     end
% 
%     printing(string)
%     close all
    
%Instantaneous Frequency
% c = categorical(cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)); 
% bar(c,fi_cortex)
% ylabel('Average instantaneous frequency')
% title(xx{1})
% ylim([100 300])
% 
%     if size(label1,1)~=3  % IF not Plusmaze 
%       string=strcat('HFOs_counts_',xx{1},'_Rat',num2str(Rat),'_',labelconditions{iii}); 
%     else
% %         if strcmp(xx{1},'HPC')
% %                   string=strcat('HFOs_counts_',xx{1},'_Rat',num2str(Rat),'_',num2str(tr(1)));         
% %         else
%                   string=strcat('HFOs_instantaneous_frequency_',xx{1},'_Rat',num2str(Rat),'_',num2str(tr(2)));         
% %         end
%     end
% 
%     printing(string)
%     close all
    
%Amplitude    
% c = categorical(cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)); 
% bar(c,amp_cortex)
% ylabel('Amplitude (uV)')
% title(xx{1})
% %ylim([100 300])
% 
%     if size(label1,1)~=3  % IF not Plusmaze 
%       string=strcat('HFOs_counts_',xx{1},'_Rat',num2str(Rat),'_',labelconditions{iii}); 
%     else
% %         if strcmp(xx{1},'HPC')
% %                   string=strcat('HFOs_counts_',xx{1},'_Rat',num2str(Rat),'_',num2str(tr(1)));         
% %         else
%                   string=strcat('HFOs_amplitude_',xx{1},'_Rat',num2str(Rat),'_',num2str(tr(2)));         
% %         end
%     end
% 
%     printing(string)
%     close all
    
    
%     xo
    TT=table;
    TT.Variables=    [[{'Count'};{'Rate'};{'Duration'};{'Average Frequency'};{'Instantaneous Frequency'};{'Amplitude'}] num2cell([hfos_cortex;hfos_cortex_rate;hfos_cortex_duration;fa_cortex;fi_cortex; amp_cortex])];
%     TT.Variables=    [[{'Count'};{'Rate'};{'Duration'}] num2cell([hfos_hpc;hfos_hpc_rate;hfos_hpc_duration])];
    
    TT.Properties.VariableNames=['Metric';cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)].';
    
%     if strcmp(xx{1},'HPC')
%             writetable(TT,strcat(xx{1},'_',num2str(tr(1)),'.xls'),'Sheet',1,'Range','A2:L6')    
%     else
            writetable(TT,strcat(xx{1},'_',num2str(tr(2)),'.xls'),'Sheet',1,'Range','A2:L10')    
%     end
%%
%Multiplets
t1=repmat({'x'},[1 length(g)+2]);

for ll=1:3
    TT=table;
%     strcat('TT.Variables=    [[','{' ,'''','Count','''','};','{' ,'''','Rate','''','};','{' ,'''','Duration','''','};',']')
    eval(strcat('TT.Variables=    [','[','{' ,'''',multiplets{ll},'''','};','{' ,'''','x','''','};','{' ,'''','x','''','}',']'," ",'[','{' ,'''','Count','''','};','{' ,'''','Rate','''','};','{' ,'''','Duration','''','}',']',...
    ' num2cell([hfos_cortex_',multiplets{ll},';hfos_cortex_rate_',multiplets{ll},';hfos_cortex_duration_',multiplets{ll},'])];'))
    TT.Properties.VariableNames=['Event';'Metric';cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)].';

    tab.(multiplets{ll})=TT;
%   eval(strcat('writetable(TT,strcat(''','HPC','''',',','''','_','''',',','num2str(tr(1)),''_',multiplets{ll},'''',',','''','.xls','''','),',...
%       '''','Sheet','''',',1,','''','Range','''',',','''','A2:L10','''',')'))
    if ll==1
        Tab=tab.(multiplets{ll});
    else
        Tab=[Tab;t1;tab.(multiplets{ll})];
    end

end

writetable(Tab,strcat(xx{1},'_',num2str(tr(2)),'_multiplets','.xls'),'Sheet',1,'Range','A1:Z50')

%%
%Cortex cohfos
    TT=table;
    TT.Variables=    [[{'Count'};{'Rate'};{'Duration'};{'Average Frequency'};{'Instantaneous Frequency'};{'Amplitude'}] num2cell([count_cohfo_cortex;rate_cohfo_cortex;dura_cohfo_cortex;fa_cohfo_cortex;fi_cohfo_cortex; amp_cohfo_cortex])];
%     TT.Variables=    [[{'Count'};{'Rate'};{'Duration'}] num2cell([hfos_hpc;hfos_hpc_rate;hfos_hpc_duration])];
    
    TT.Properties.VariableNames=['Metric';cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)].';
    
%     if strcmp(xx{1},'HPC')
%             writetable(TT,strcat(xx{1},'_',num2str(tr(1)),'.xls'),'Sheet',1,'Range','A2:L6')    
%     else
            writetable(TT,strcat(xx{1},'_',num2str(tr(2)),'_cohfos','.xls'),'Sheet',1,'Range','A2:L10')    

%Cortex singles
    TT=table;
    TT.Variables=    [[{'Count'};{'Rate'};{'Duration'};{'Average Frequency'};{'Instantaneous Frequency'};{'Amplitude'}] num2cell([count_single_cortex;rate_single_cortex;dura_single_cortex;fa_single_cortex;fi_single_cortex; amp_single_cortex])];
%     TT.Variables=    [[{'Count'};{'Rate'};{'Duration'}] num2cell([hfos_hpc;hfos_hpc_rate;hfos_hpc_duration])];
    
    TT.Properties.VariableNames=['Metric';cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)].';
    
%     if strcmp(xx{1},'HPC')
%             writetable(TT,strcat(xx{1},'_',num2str(tr(1)),'.xls'),'Sheet',1,'Range','A2:L6')    
%     else
            writetable(TT,strcat(xx{1},'_',num2str(tr(2)),'_singles','.xls'),'Sheet',1,'Range','A2:L10')    

  %%          

%HPC
c = categorical(cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)); 
bar(c,hfos_hpc)
ylabel('Number of HFOs')
title('HPC')

    if size(label1,1)~=3  % IF not Plusmaze 
      string=strcat('HFOs_counts_',xx{1},'_Rat',num2str(Rat),'_',labelconditions{iii}); 
    else
%         if strcmp(xx{1},'HPC')
%                   string=strcat('HFOs_counts_',xx{1},'_Rat',num2str(Rat),'_',num2str(tr(1)));         
%         else
                  string=strcat('HFOs_counts_','HPC','_Rat',num2str(Rat),'_',num2str(tr(1)));         
%         end
    end

    printing(string)
    close all
%rate
c = categorical(cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)); 
bar(c,hfos_hpc_rate)
ylabel('HFOs per second')
title('HPC')


    if size(label1,1)~=3  % IF not Plusmaze 
      string=strcat('HFOs_rate_',xx{1},'_Rat',num2str(Rat),'_',labelconditions{iii}); 
    else
%         if strcmp(xx{1},'HPC')
%                   string=strcat('HFOs_rate_',xx{1},'_Rat',num2str(Rat),'_',num2str(tr(1)));         
%         else
                  string=strcat('HFOs_rate_','HPC','_Rat',num2str(Rat),'_',num2str(tr(1)));         
%         end
        
    end

    printing(string)
    close all
%    xo


%Average Frequency
c = categorical(cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)); 
bar(c,fa_hpc)
ylabel('Average frequency')
title('HPC')
ylim([100 300])

    if size(label1,1)~=3  % IF not Plusmaze 
      string=strcat('HFOs_counts_',xx{1},'_Rat',num2str(Rat),'_',labelconditions{iii}); 
    else
%         if strcmp(xx{1},'HPC')
%                   string=strcat('HFOs_counts_',xx{1},'_Rat',num2str(Rat),'_',num2str(tr(1)));         
%         else
                  string=strcat('HFOs_average_frequency_','HPC','_Rat',num2str(Rat),'_',num2str(tr(1)));         
%         end
    end

    printing(string)
    close all
    
%Instantaneous Frequency
c = categorical(cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)); 
bar(c,fi_hpc)
ylabel('Average instantaneous frequency')
title('HPC')
ylim([100 300])

    if size(label1,1)~=3  % IF not Plusmaze 
      string=strcat('HFOs_counts_',xx{1},'_Rat',num2str(Rat),'_',labelconditions{iii}); 
    else
%         if strcmp(xx{1},'HPC')
%                   string=strcat('HFOs_counts_',xx{1},'_Rat',num2str(Rat),'_',num2str(tr(1)));         
%         else
                  string=strcat('HFOs_instantaneous_frequency_','HPC','_Rat',num2str(Rat),'_',num2str(tr(1)));         
%         end
    end

    printing(string)
    close all

%Amplitude    
c = categorical(cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)); 
bar(c,amp_hpc)
ylabel('Amplitude (uV)')
title('HPC')
%ylim([100 300])

    if size(label1,1)~=3  % IF not Plusmaze 
      string=strcat('HFOs_counts_',xx{1},'_Rat',num2str(Rat),'_',labelconditions{iii}); 
    else
%         if strcmp(xx{1},'HPC')
%                   string=strcat('HFOs_counts_',xx{1},'_Rat',num2str(Rat),'_',num2str(tr(1)));         
%         else
                  string=strcat('HFOs_amplitude_','HPC','_Rat',num2str(Rat),'_',num2str(tr(1)));         
%         end
    end

    printing(string)
    close all


    TT=table;
    TT.Variables=    [[{'Count'};{'Rate'};{'Duration'};{'Average Frequency'};{'Instantaneous Frequency'};{'Amplitude'}] num2cell([hfos_hpc;hfos_hpc_rate;hfos_hpc_duration;fa_hpc;fi_hpc;amp_hpc])];        
%    TT.Properties.VariableNames=['Metric';g];
    TT.Properties.VariableNames=['Metric';cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)].';
    
%     if strcmp(xx{1},'HPC')
%             writetable(TT,strcat(xx{1},'_',num2str(tr(1)),'.xls'),'Sheet',1,'Range','A2:L6')    
%     else
            writetable(TT,strcat('HPC','_',num2str(tr(1)),'.xls'),'Sheet',1,'Range','A2:L10')    
%     end
%%
%Multiplets
t1=repmat({'x'},[1 length(g)+2]);

for ll=1:length(multiplets)
    TT=table;
%     strcat('TT.Variables=    [[','{' ,'''','Count','''','};','{' ,'''','Rate','''','};','{' ,'''','Duration','''','};',']')
    eval(strcat('TT.Variables=    [','[','{' ,'''',multiplets{ll},'''','};','{' ,'''','x','''','};','{' ,'''','x','''','}',']'," ",'[','{' ,'''','Count','''','};','{' ,'''','Rate','''','};','{' ,'''','Duration','''','}',']',...
    ' num2cell([hfos_hpc_',multiplets{ll},';hfos_hpc_rate_',multiplets{ll},';hfos_hpc_duration_',multiplets{ll},'])];'))
    TT.Properties.VariableNames=['Event';'Metric';cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)].';

    tab.(multiplets{ll})=TT;
%   eval(strcat('writetable(TT,strcat(''','HPC','''',',','''','_','''',',','num2str(tr(1)),''_',multiplets{ll},'''',',','''','.xls','''','),',...
%       '''','Sheet','''',',1,','''','Range','''',',','''','A2:L10','''',')'))
    if ll==1
        Tab=tab.(multiplets{ll});
    else
        Tab=[Tab;t1;tab.(multiplets{ll})];
    end

end

writetable(Tab,strcat('HPC','_',num2str(tr(1)),'_multiplets','.xls'),'Sheet',1,'Range','A1:Z50')


%%
%hpc cohfos
    TT=table;
    TT.Variables=    [[{'Count'};{'Rate'};{'Duration'};{'Average Frequency'};{'Instantaneous Frequency'};{'Amplitude'}] num2cell([count_cohfo_hpc;rate_cohfo_hpc;dura_cohfo_hpc;fa_cohfo_hpc;fi_cohfo_hpc; amp_cohfo_hpc])];
%     TT.Variables=    [[{'Count'};{'Rate'};{'Duration'}] num2cell([hfos_hpc;hfos_hpc_rate;hfos_hpc_duration])];
    
    TT.Properties.VariableNames=['Metric';cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)].';
    
%     if strcmp(xx{1},'HPC')
%             writetable(TT,strcat(xx{1},'_',num2str(tr(1)),'.xls'),'Sheet',1,'Range','A2:L6')    
%     else
            writetable(TT,strcat('HPC','_',num2str(tr(1)),'_',num2str(tr(2)),'_cohfos','.xls'),'Sheet',1,'Range','A2:L10')    

%hpc singles
    TT=table;
    TT.Variables=    [[{'Count'};{'Rate'};{'Duration'};{'Average Frequency'};{'Instantaneous Frequency'};{'Amplitude'}] num2cell([count_single_hpc;rate_single_hpc;dura_single_hpc;fa_single_hpc;fi_single_hpc; amp_single_hpc])];
%     TT.Variables=    [[{'Count'};{'Rate'};{'Duration'}] num2cell([hfos_hpc;hfos_hpc_rate;hfos_hpc_duration])];
    
    TT.Properties.VariableNames=['Metric';cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)].';
    
%     if strcmp(xx{1},'HPC')
%             writetable(TT,strcat(xx{1},'_',num2str(tr(1)),'.xls'),'Sheet',1,'Range','A2:L6')    
%     else
            writetable(TT,strcat('HPC','_',num2str(tr(1)),'_',num2str(tr(2)),'_singles','.xls'),'Sheet',1,'Range','A2:L10')    



%%
%COHFOS
c = categorical(cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)); 
bar(c,cohfos_count)
ylabel('Number of coHFOs')
title('Both areas')

    if size(label1,1)~=3  % IF not Plusmaze 
      string=strcat('coHFOs_counts_',xx{1},'_Rat',num2str(Rat),'_',labelconditions{iii}); 
    else
%         if strcmp(xx{1},'HPC')
%                   string=strcat('HFOs_counts_',xx{1},'_Rat',num2str(Rat),'_',num2str(tr(1)));         
%         else
                  string=strcat('coHFOs_counts_','Rat',num2str(Rat),'_',num2str(tr(1)),'_',num2str(tr(2)));         
%         end
    end

    printing(string)
    close all
    
    
c = categorical(cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)); 
bar(c,cohfos_rate)
ylabel('coHFOs per second')
title('Both areas')    

    if size(label1,1)~=3  % IF not Plusmaze 
      string=strcat('coHFOs_counts_',xx{1},'_Rat',num2str(Rat),'_',labelconditions{iii}); 
    else
%         if strcmp(xx{1},'HPC')
%                   string=strcat('HFOs_counts_',xx{1},'_Rat',num2str(Rat),'_',num2str(tr(1)));         
%         else
                  string=strcat('coHFOs_rate_','Rat',num2str(Rat),'_',num2str(tr(1)),'_',num2str(tr(2)));         
%         end
    end

    printing(string)
    close all

    TT=table;
    TT.Variables=    [[{'Count'};{'Rate'}] num2cell([cohfos_count;cohfos_rate;])];
    TT.Properties.VariableNames=['Metric';g];    
    writetable(TT,strcat('coHFOs_',num2str(tr(2)),'.xls'),'Sheet',1,'Range','A2:L6')    
    
%%
%Multiplet cohfos
t1=repmat({'x'},[1 length(g)+2]);
for ll=1:length(multiplets)

    TT=table;
    TT.Variables=    [[{multiplets{ll}};{'x'}] [{'Count'};{'Rate'}] num2cell([cohfos_count_multiplets.(multiplets{ll});cohfos_rate_multiplets.(multiplets{ll});])];
  
    TT.Properties.VariableNames=['Event';'Metric';g];
    tab_cohfos.(multiplets{ll})=TT;
        if ll==1
            Tab_cohfos=tab_cohfos.(multiplets{ll});
        else
            Tab_cohfos=[Tab_cohfos;t1;tab_cohfos.(multiplets{ll})];
        end

%     writetable(TT,strcat('coHFOs_',num2str(tr(2)),'.xls'),'Sheet',1,'Range','A2:L6')
end


writetable(Tab_cohfos,strcat('coHFOs','_multiplets_',num2str(tr(1)),'_',num2str(tr(2)),'.xls'),'Sheet',1,'Range','A1:Z50')

%% Demixed Gaussians 

%coHFOs
    TT=table;
%    TT.Variables=    [[{'Slower'};{'x'}] [{'Count'};{'Rate'}] num2cell([cohfos_count_g1;cohfos_rate_g1;])];
    TT.Variables=    [[{'Slower'};{'x'};{'Faster'};{'x'}] [{'Count'};{'Rate'};{'Count'};{'Rate'}] num2cell([cohfos_count_g1;cohfos_rate_g1;cohfos_count_g2;cohfos_rate_g2])];

    TT.Properties.VariableNames=['Events';'Metric';g];


writetable(TT,strcat('SWRD4_slower_faster_cohfos_',num2str(tr(1)),'_',num2str(tr(2)),'.xls'),'Sheet',1,'Range','A1:Z50')

%Singles
    TT=table;
%    TT.Variables=    [[{'Slower'};{'x'}] [{'Count'};{'Rate'}] num2cell([cohfos_count_g1;cohfos_rate_g1;])];
    TT.Variables=    [[{'Slower'};{'x'};{'Faster'};{'x'}] [{'Count'};{'Rate'};{'Count'};{'Rate'}] num2cell([singles_count_g1;singles_rate_g1;singles_count_g2;singles_rate_g2])];

    TT.Properties.VariableNames=['Events';'Metric';g];


writetable(TT,strcat('SWRD4_slower_faster_singles_',num2str(tr(1)),'_',num2str(tr(2)),'.xls'),'Sheet',1,'Range','A1:Z50')

%% Sum total of events

%Counts
    TT=table;
%     TT.Variables=    [[{'Slower'};{'x'};{'Faster'};{'x'}] [{'Count'};{'Rate'};{'Count'};{'Rate'}] num2cell([singles_count_g1;singles_rate_g1;singles_count_g2;singles_rate_g2])];
    TT.Variables=    [[{'Cohfos slow'};{'Cohfos fast'};{'Single slow'};{'Single fast'};{'All slow'};{'All fast'};{'Total'}] num2cell([cohfos_count_g1;cohfos_count_g2;singles_count_g1;singles_count_g2;all_count_g1;all_count_g2;all_count])];

    TT.Properties.VariableNames=['Events';g];

    writetable(TT,strcat('SWRD4_counts_',num2str(tr(1)),'_',num2str(tr(2)),'.xls'),'Sheet',1,'Range','A1:Z50')

    
%Rates
    TT=table;
    TT.Variables=    [[{'Cohfos slow'};{'Cohfos fast'};{'Single slow'};{'Single fast'};{'All slow'};{'All fast'};{'Total'}] num2cell([cohfos_rate_g1;cohfos_rate_g2;singles_rate_g1;singles_rate_g2;all_rate_g1;all_rate_g2;all_rate])];
% 
%     TT.Variables=    [[{'All slow'};{'All fast'};{'Total'}] num2cell([all_rate_g1;all_rate_g2;all_rate])];

    TT.Properties.VariableNames=['Events';g];
    writetable(TT,strcat('SWRD4_rates_',num2str(tr(1)),'_',num2str(tr(2)),'.xls'),'Sheet',1,'Range','A1:Z50')

    %%
    %%
all_count_g1
all_rate_g1

all_count_g2
all_rate_g2

all_count
all_rate
xo