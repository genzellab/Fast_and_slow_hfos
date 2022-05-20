% GL_delta_count.m
%Detects Delta waves and counts them.
% Requires 'load_me_first.mat' loaded first. 
clear variables
cd('/home/adrian/Documents/GitHub/CorticoHippocampal/Fast_and_slow_hfos')
load('load_me_first.mat')

%% Find location
close all

% dname=uigetdir([],'Select folder with Matlab data containing all rats.');
% cd(dname)

%cd('/home/adrian/Documents/Plusmaze_downsampled')

dname='/media/adrian/6aa1794c-0320-4096-a7df-00ab0ba946dc/Plusmaze_downsampled/Data_plusmaze';
cd(dname)

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
%     error('Name issue')
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

pfc_thresholds=[1.5, 3, -1.5, 0];

%Find PFC Deltawaves

[deltaWave_count_pfc,deltaFreq_pfc,delta_duration_pfc,Mx_pfc,timeasleep,sig_pfc,Ex_pfc,Sx_pfc, ...
   DeltaWaves, ti_cont,duration_epoch_cumsum]=gui_finddeltawavesZugaro(CORTEX,states,xx,multiplets,fn, pfc_thresholds);

deltaFreq_PFC(k)=deltaFreq_pfc;
delta_duration_PFC(k)=delta_duration_pfc;
delta_count_PFC(k)=deltaWave_count_pfc;

sig_pfc =[sig_pfc{:}];
print_hist=0;

%This delta_specs comes from Milan's RGS repo.
[x,y,z,w,h,q,l,p,si_mixed,th,PCA_features]=delta_specs(sig_pfc,timeasleep,print_hist);

instafreq_pfc{k}=x;
avgfreq_pfc{k}=y;
magnitude_pfc{k}=z;
auc_pfc{k}=l;
duration_pfc{k}=q;
p2p_pfc{k}=p;


%Find PPC Deltawaves
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

par_thresholds= [1.5, 3, -1.5, 0];

[deltaWave_count_par,deltaFreq_par,delta_duration_par,Mx_par,timeasleep,sig_par,Ex_par,Sx_par, ...
   DeltaWaves, ti_cont,duration_epoch_cumsum]=gui_finddeltawavesZugaro(par,states,yy,multiplets,fn, par_thresholds);

deltaFreq_PAR(k)=deltaFreq_par;
delta_duration_PAR(k)=delta_duration_par;
delta_count_PAR(k)=deltaWave_count_par;

sig_par =[sig_par{:}];

%This delta_specs comes from Milan's RGS repo.
[x,y,z,w,h,q,l,p,si_mixed,th,PCA_features]=delta_specs(sig_par,timeasleep,print_hist);

instafreq_par{k}=x;
avgfreq_par{k}=y;
magnitude_par{k}=z;
auc_par{k}=l;
duration_par{k}=q;
p2p_par{k}=p;

progress_bar(k,length(g),f)
cd ..    
end

% plotting all delta waves
figure,
for i=1:deltaWave_count_pfc
    plot(sig_pfc{i})
    hold on;
end
xlabel('duration')
ylabel('amplitute')

% plotting PFC specs
figure,
for k=1:4
    i=5-k;
    subplot(3,2,1)
    hold on;
    histogram(instafreq_pfc{i},'BinWidth', 0.5); 
    title('Instantaneous Frequencies');xlabel('Frequency (Hz)');ylabel('Count') 
    hold off;
    legend(cell2mat(labelconditions2(4)), cell2mat(labelconditions2(3)), cell2mat(labelconditions2(2)), cell2mat(labelconditions2(1)) )
    
    subplot(3,2,2)
    hold on;
    histogram(avgfreq_pfc{i}, 'BinWidth', 0.25);
    title('Average Frequencies');xlabel('Frequency (Hz)');ylabel('Count') 
    hold off;
    legend(cell2mat(labelconditions2(4)), cell2mat(labelconditions2(3)), cell2mat(labelconditions2(2)), cell2mat(labelconditions2(1)) )
        
    subplot(3,2,3)
    hold on;
    histogram(magnitude_pfc{i}, 'BinWidth', 5);
    title('Amplitude');xlabel('\muV');ylabel('Count')   
    hold off;
    legend(cell2mat(labelconditions2(4)), cell2mat(labelconditions2(3)), cell2mat(labelconditions2(2)), cell2mat(labelconditions2(1)) )
        
    subplot(3,2,4)
    hold on;
    histogram(auc_pfc{i}, 'BinWidth', 0.5);
    title('Area under the curve');xlabel('AUC');ylabel('Count') 
    hold off;
    legend(cell2mat(labelconditions2(4)), cell2mat(labelconditions2(3)), cell2mat(labelconditions2(2)), cell2mat(labelconditions2(1)) )
        
    subplot(3,2,5)
    hold on;
    histogram(duration_pfc{i},  'BinWidth', 0.01);
    title('Duration');xlabel('Miliseconds');ylabel('Count')   
    hold off;
    legend(cell2mat(labelconditions2(4)), cell2mat(labelconditions2(3)), cell2mat(labelconditions2(2)), cell2mat(labelconditions2(1)) )
    
    subplot(3,2,6)
    hold on;
    histogram(p2p_pfc{i}, 'BinWidth', 5);
    title('Peak-to-peak amplitude');xlabel('\muV');ylabel('Count');   
    hold off;
    legend(cell2mat(labelconditions2(4)), cell2mat(labelconditions2(3)), cell2mat(labelconditions2(2)), cell2mat(labelconditions2(1)) )
        
end

% plotting PAR specs
figure,
for k=1:4
    i=5-k;
    subplot(3,2,1)
    hold on;
    histogram(instafreq_par{i}, 'BinWidth', 0.5); 
    title('Instantaneous Frequencies');xlabel('Frequency (Hz)');ylabel('Count') 
    hold off;
    legend(cell2mat(labelconditions2(4)), cell2mat(labelconditions2(3)), cell2mat(labelconditions2(2)), cell2mat(labelconditions2(1)) )
    
    subplot(3,2,2)
    hold on;
    histogram(avgfreq_par{i}, 'BinWidth', 0.25);
    title('Average Frequencies');xlabel('Frequency (Hz)');ylabel('Count') 
    hold off;
    legend(cell2mat(labelconditions2(4)), cell2mat(labelconditions2(3)), cell2mat(labelconditions2(2)), cell2mat(labelconditions2(1)) )
        
    subplot(3,2,3)
    hold on;
    histogram(magnitude_par{i},  'BinWidth', 5);
    title('Amplitude');xlabel('\muV');ylabel('Count')   
    hold off;
    legend(cell2mat(labelconditions2(4)), cell2mat(labelconditions2(3)), cell2mat(labelconditions2(2)), cell2mat(labelconditions2(1)) )
        
    subplot(3,2,4)
    hold on;
    histogram(auc_par{i}, 'BinWidth', 0.5);
    title('Area under the curve');xlabel('AUC');ylabel('Count') 
    hold off;
    legend(cell2mat(labelconditions2(4)), cell2mat(labelconditions2(3)), cell2mat(labelconditions2(2)), cell2mat(labelconditions2(1)) )
        
    subplot(3,2,5)
    hold on;
    histogram(duration_par{i},  'BinWidth', 0.01);
    title('Duration');xlabel('Miliseconds');ylabel('Count')   
    hold off;
    legend(cell2mat(labelconditions2(4)), cell2mat(labelconditions2(3)), cell2mat(labelconditions2(2)), cell2mat(labelconditions2(1)) )
    
    subplot(3,2,6)
    hold on;
    histogram(p2p_par{i},  'BinWidth', 5);
    title('Peak-to-peak amplitude');xlabel('\muV');ylabel('Count');   
    hold off;
    legend(cell2mat(labelconditions2(4)), cell2mat(labelconditions2(3)), cell2mat(labelconditions2(2)), cell2mat(labelconditions2(1)) )
        
end

xo
cd ..

% table
excel_range=struct('rat_24', {'A2:L7' 'A1'}, 'rat_26', {'A8:L16' 'A7'}, 'rat_27', {'A14:L22' 'A13'});
rat_no=strcat('rat_', num2str(Rat));

TT=table;
TT.Variables= [[{'Count'};{'Rate'};{'Duration'}] num2cell([delta_count_PFC; deltaFreq_PFC; delta_duration_PFC])];
TT.Properties.VariableNames=['Metric' cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)].';
writematrix('PFC','DeltaWave.xls','Sheet',1,'Range',excel_range(2).(rat_no))   
writetable(TT,'DeltaWave.xls','Sheet',1,'Range',excel_range(1).(rat_no))    

TT=table;
TT.Variables= [[{'Count'};{'Rate'};{'Duration'}] num2cell([delta_count_PAR; deltaFreq_PAR; delta_duration_PAR])];
TT.Properties.VariableNames=['Metric' cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)].';
writematrix('PAR','DeltaWave.xls','Sheet',2, 'Range', excel_range(2).(rat_no))   
writetable(TT,'DeltaWave.xls','Sheet',2, 'Range', excel_range(1).(rat_no))    


