%GL_delta_spindle
%Detects spindles and computes coocurrence with delta waves.
% Requires 'load_me_first.mat' loaded first. 
clear variables
cd('/home/adrian/Documents/GitHub/CorticoHippocampal/Fast_and_slow_hfos')
load('load_me_first.mat')

%% Find location
close all

% dname=uigetdir([],'Select folder with Matlab data containing all rats.');
% cd(dname)

% cd('/home/adrian/Documents/Plusmaze_downsampled')

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
% if sum(cell2mat(cellfun(@(equis1) contains(equis1,'nl'),g,'UniformOutput',false)))==1
% g=g([find(cell2mat(cellfun(@(equis1) contains(equis1,labelconditions2{1}),g,'UniformOutput',false)))...
%  find(cell2mat(cellfun(@(equis1) contains(equis1,labelconditions2{2}),g,'UniformOutput',false)))...
%  find(cell2mat(cellfun(@(equis1) contains(equis1,labelconditions2{3}),g,'UniformOutput',false)))...
%  find(cell2mat(cellfun(@(equis1) contains(equis1,labelconditions2{4}),g,'UniformOutput',false)))]);
% 
% else
%     error('Name issue')
% end

%Get thresholds for event detection.
tr=getfield(T,strcat('Rat',num2str(Rat)));%Thresholds 

tic
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
%xo
% PFC spindles
[spindle_count_pfc,~,~,Mx_cortex,timeasleep,sig_cortex,Ex_cortex,Sx_cortex]=gui_findspindlesZugaro(CORTEX,states,xx,multiplets,fn);
%xo


%Find PFC Deltawaves
pfc_thresholds=[1.5, 3, -1.5, 0];
[deltaWave_count_pfc,deltaFreq_pfc,delta_duration_pfc,Mx_pfc,timeasleep,sig_pfc,Ex_pfc,Sx_pfc, ...
   DeltaWaves, ti_cont,duration_epoch_cumsum]=gui_finddeltawavesZugaro(CORTEX,states,xx,multiplets,fn, pfc_thresholds);


%% PAR Spindles and deltawaves

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
[spindle_count_par,RipFreq,rip_duration,Mx_par,timeasleep,sig_par,Ex_par,Sx_par,...
    ]=gui_findspindlesZugaro(par,states,yy,multiplets,fn);

%PAR deltawaves
par_thresholds= [1.5, 3, -1.5, 0];
[deltaWave_count_par,deltaFreq_par,delta_duration_par,Mx_par_delta,timeasleep,sig_par_delta,Ex_par_delta,Sx_par_delta, ...
   DeltaWaves, ti_cont,duration_epoch_cumsum]=gui_finddeltawavesZugaro(par,states,yy,multiplets,fn, par_thresholds);


%Ex_cortex: PFC SPINDLE
%Ex_pfc:   PFC Deltawaves
%Ex_par:  PAR SPINDLE
%Ex_par_delta: PAR Deltawaves

%% Coocur PFC spindle and PFC delta
[cohfos1_deltaPFC_spindlePFC,cohfos2_deltaPFC_spindlePFC]=cellfun(@(equis1,equis2) co_hfo_delta_spindle(equis1,equis2),Mx_pfc.',Mx_cortex,'UniformOutput',false);
Cooccur_deltaPFC_spindlePFC{k}=([cohfos1_deltaPFC_spindlePFC{:}]);
Count_deltaPFC_spindle_PFC(k)=length([cohfos1_deltaPFC_spindlePFC{:}]);
Rate_deltaPFC_spindle_PFC(k)=length([cohfos1_deltaPFC_spindlePFC{:}])/timeasleep;

%% Coocur PFC spindle and PAR delta
[cohfos1_deltaPAR_spindlePFC,cohfos2_deltaPAR_spindlePFC]=cellfun(@(equis1,equis2) co_hfo_delta_spindle(equis1,equis2),Mx_par_delta.',Mx_cortex,'UniformOutput',false);
Cooccur_deltaPAR_spindlePFC{k}=([cohfos1_deltaPAR_spindlePFC{:}]);
Count_deltaPAR_spindlePFC(k)=length([cohfos1_deltaPAR_spindlePFC{:}]);
Rate_deltaPAR_spindlePFC(k)=length([cohfos1_deltaPAR_spindlePFC{:}])/timeasleep;

%% Coocur PAR spindle and PFC delta
[cohfos1_deltaPFC_spindlePAR,cohfos2_deltaPFC_spindlePAR]=cellfun(@(equis1,equis2) co_hfo_delta_spindle(equis1,equis2),Mx_pfc.',Mx_par,'UniformOutput',false);
Cooccur_deltaPFC_spindlePAR{k}=([cohfos1_deltaPFC_spindlePAR{:}]);
Count_deltaPFC_spindlePAR(k)=length([cohfos1_deltaPFC_spindlePAR{:}]);
Rate_deltaPFC_spindlePAR(k)=length([cohfos1_deltaPFC_spindlePAR{:}])/timeasleep;

%% Coocur PAR spindle and PAR delta
[cohfos1_deltaPAR_spindlePAR,cohfos2_deltaPAR_spindlePAR]=cellfun(@(equis1,equis2) co_hfo_delta_spindle(equis1,equis2),Mx_par_delta.',Mx_par,'UniformOutput',false);
Cooccur_deltaPAR_spindlePAR{k}=([cohfos1_deltaPAR_spindlePAR{:}]);
Count_deltaPAR_spindlePAR(k)=length([cohfos1_deltaPAR_spindlePAR{:}]);
Rate_deltaPAR_spindlePAR(k)=length([cohfos1_deltaPAR_spindlePAR{:}])/timeasleep;


% PAR Spindle - PAR delta
% PAR Spindle - PFC delta
% PFC Spindle - PAR delta
% PFC Spindle - PFC delta



progress_bar(k,length(g),f)
    cd ..    
    end
toc
%%
%xo
Counts=[Count_deltaPFC_spindle_PFC;...
Count_deltaPAR_spindlePFC;...
Count_deltaPFC_spindlePAR;...
Count_deltaPAR_spindlePAR]

Rates=[Rate_deltaPFC_spindle_PFC;...
Rate_deltaPAR_spindlePFC;...
Rate_deltaPFC_spindlePAR;...
Rate_deltaPAR_spindlePAR]

TT=table;
TT.Variables=[Counts;Rates];
TT.Properties.VariableNames=[{'baseline'},{'foraging'},{'novelty'},{'plusmaze'}]
writetable(TT,'DELTA_SPINDLE.xls','Sheet',1,'Range','A2:Z50')
% save('spindle-HFO-ripple Coocurances.mat', 'Cohfos1_PAR_g1_pre','Cohfos1_PAR_g2_pre','Cohfos1_PFC_g1_pre','Cohfos1_PFC_g2_pre',...
%     'Cohfos1_PAR_g1_post','Cohfos1_PAR_g2_post','Cohfos1_PFC_g1_post','Cohfos1_PFC_g2_post',...
%     'Cohfos2_PAR_hpc_pre','Cohfos2_PAR_hpc_post','Cohfos2_PFC_hpc_pre','Cohfos2_PFC_hpc_post',...
%     'Out_PAR_pre','Out_PFC_pre','Out_PAR_post','Out_PFC_post');
xo
cd ..        
%% Generate tables and save values into spreadsheets.

excel_range=struct('rat_24', {'A2:L10' 'A1' 'A2:L39' 'A1'}, 'rat_26', {'A11:L18' 'A10' 'A40:L77' 'A39'}, 'rat_27', {'A20:L26' 'A19' 'A77:L114' 'A76'});
rat_no=strcat('rat_', num2str(Rat));

  %% HFOs-spindles PRE and POST
    TT=table;
    TT.Variables= [[{'PAR slow pre'};{'PAR fast Pre'};{'PFC slow pre'};{'PFC Fast pre'}], num2cell([cellfun('length',Cohfos1_PAR_g1_pre); cellfun('length',Cohfos1_PAR_g2_pre);...
       cellfun('length',Cohfos1_PFC_g1_pre); cellfun('length',Cohfos1_PFC_g2_pre)])];
   
    TT.Properties.VariableNames=['Metric' cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)];
    
    writematrix(rat_no,strcat('ripple_spindle_HFO','.xls'),'Sheet',1,'Range',excel_range(2).(rat_no))
    writematrix('HFOs-spindles PRE',strcat('ripple_spindle_HFO','.xls'), 'sheet',1, 'Range', 'B1')
    writetable(TT,strcat('ripple_spindle_HFO','.xls'), 'Sheet',1,'Range',excel_range(1).(rat_no))    

    TT=table;
    TT.Variables= [[{'PAR slow post'};{'PAR fast Post'};{'PFC slow post'};{'PFC Fast post'}] num2cell([cellfun('length',Cohfos1_PAR_g1_post); cellfun('length',Cohfos1_PAR_g2_post);...
       cellfun('length',Cohfos1_PFC_g1_post); cellfun('length',Cohfos1_PFC_g2_post)])]; 

    TT.Properties.VariableNames=['Metric' cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)];
    
    writematrix(rat_no,strcat('ripple_spindle_HFO','.xls'),'Sheet',2,'Range',excel_range(2).(rat_no))
    writematrix('HFOs-spindles POST',strcat('ripple_spindle_HFO','.xls'), 'sheet',2, 'Range', 'B1')
    writetable(TT,strcat('ripple_spindle_HFO','.xls'),'Sheet',2,'Range',excel_range(1).(rat_no))    

%% HPC ripples-spindles PRE POST
    TT=table;
    TT.Variables= [[{'PAR all pre'};{'PAR Unique Pre'};{'PFC all pre'};{'PFC unique pre'}], num2cell([cellfun('length',Cohfos2_PAR_hpc_pre); cell2mat(cellfun(@(equis) length(unique(equis)),Cohfos2_PAR_hpc_pre,'UniformOutput',false));...
       cellfun('length',Cohfos2_PFC_hpc_pre); cell2mat(cellfun(@(equis) length(unique(equis)),Cohfos2_PFC_hpc_pre,'UniformOutput',false)) ])]; 
    
    TT.Properties.VariableNames=['Metric' cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)];
    
                writematrix(rat_no,strcat('ripple_spindle_HFO','.xls'),'Sheet',3,'Range',excel_range(2).(rat_no))
                writematrix('HPC Ripples-spindles PRE', strcat('ripple_spindle_HFO','.xls'),'sheet',3, 'Range', 'B1')
                writetable(TT,strcat('ripple_spindle_HFO','.xls'),'Sheet',3,'Range',excel_range(1).(rat_no))    

                    TT=table;
    TT.Variables= [[{'PAR all Post'};{'PAR Unique Post'};{'PFC all post'};{'PFC unique post'}] num2cell([cellfun('length',Cohfos2_PAR_hpc_post); ...
                    cell2mat(cellfun(@(equis) length(unique(equis)), Cohfos2_PAR_hpc_post,'UniformOutput',false));...
                    cellfun('length',Cohfos2_PFC_hpc_post); ...
                    cell2mat(cellfun(@(equis) length(unique(equis)),Cohfos2_PFC_hpc_post,'UniformOutput',false)) ]) ]; 
    
    TT.Properties.VariableNames=['Metric' cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)];
    
                writematrix(rat_no,strcat('ripple_spindle_HFO','.xls'),'Sheet',4,'Range',excel_range(2).(rat_no))
                writematrix('HPC Ripples-spindles POST',strcat('ripple_spindle_HFO','.xls'), 'sheet',4, 'Range', 'B1')
                writetable(TT,strcat('ripple_spindle_HFO','.xls'),'Sheet',4,'Range',excel_range(1).(rat_no))    

%% HPC ripples-spindles-multiplets PRE POST
    TT=table;
    TT.Variables=[[{'PAR - Single'};{'PAR - doublet'};{'PAR - triplet'};{'PAR - quadruplet '}; {'PAR- Pentuplets'};{'PAR-Sextuplets'}; {' '}; {'PFC - Single'};{'PFC - doublet'};{'PFC - triplet'};{'PFC - quadruplet '}; {'PFC-pentplets'};{'PFC-Sextuplets'}], ...
                    num2cell([ Out_PAR_pre{1}(:,2) Out_PAR_pre{2}(:,2) Out_PAR_pre{3}(:,2) Out_PAR_pre{4}(:,2); [NaN NaN NaN NaN] ;Out_PFC_pre{1}(:,2) Out_PFC_pre{2}(:,2) ...
                    Out_PFC_pre{3}(:,2) Out_PFC_pre{4}(:,2) ])];
    TT.Properties.VariableNames=['Metric' cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)];
    
                 writematrix(rat_no,strcat('ripple_spindle_HFO','.xls'),'Sheet',5,'Range',excel_range(4).(rat_no))
                 writematrix('HPC Ripples-spindles multiplets PRE',strcat('ripple_spindle_HFO','.xls'), 'sheet',5, 'Range', 'B1')
                 writetable(TT,strcat('ripple_spindle_HFO','.xls'),'Sheet',5,'Range',excel_range(3).(rat_no))    

    TT=table;
    TT.Variables= [[{'PAR - Single'};{'PAR - doublet'};{'PAR - triplet'};{'PAR - quadruplet '}; {'PAR- Pentuplets'};{'PAR-Sextuplets'}; {' '}; {'PFC - Single'};{'PFC - doublet'};{'PFC - triplet'};{'PFC - quadruplet '}; {'PFC-pentplets'};{'PFC-Sextuplets'}], ...
 num2cell([([Out_PAR_post{1}(:,2) Out_PAR_post{2}(:,2) Out_PAR_post{3}(:,2) Out_PAR_post{4}(:,2); [NaN NaN NaN NaN] ;Out_PFC_post{1}(:,2) Out_PFC_post{2}(:,2) Out_PFC_post{3}(:,2) Out_PFC_post{4}(:,2)] ) ])];
    TT.Properties.VariableNames=['Metric' cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)];
    
                 writematrix(rat_no,strcat('ripple_spindle_HFO','.xls'),'Sheet',6,'Range',excel_range(4).(rat_no))
                 writematrix('HPC Ripples-spindles multiplets POST',strcat('ripple_spindle_HFO','.xls'), 'sheet',6, 'Range', 'B1')
                 writetable(TT,strcat('ripple_spindle_HFO','.xls'),'Sheet',6,'Range',excel_range(3).(rat_no))    

  %% HFOs-spindles
    TT=table;
    TT.Variables= [[{'PAR slow'};{'PAR fast'};{'PFC slow'};{'PFC fast'}], num2cell([cellfun('length',Cohfos1_PAR_g1); cellfun('length',Cohfos1_PAR_g2);...
       cellfun('length',Cohfos1_PFC_g1); cellfun('length',Cohfos1_PFC_g2) ])]; 
   
    TT.Properties.VariableNames=['Metric' cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)];
    
            writematrix(rat_no,strcat('ripple_spindle_HFO','.xls'),'Sheet',7,'Range',excel_range(2).(rat_no))
            writematrix('HFO-spindles',strcat('ripple_spindle_HFO','.xls'), 'sheet',7, 'Range', 'B1')
            writetable(TT,strcat('ripple_spindle_HFO','.xls'),'Sheet',7,'Range', excel_range(1).(rat_no))    
%% HPC ripples-spindles
    TT=table;
    TT.Variables= [ [{'PAR all'};{'PAR unique'};{'PFC all'};{'PFC unique'}], num2cell([ cellfun('length',Cohfos2_PAR_hpc); cell2mat(cellfun(@(equis) length(unique(equis)),Cohfos2_PAR_hpc,'UniformOutput',false));...
       cellfun('length',Cohfos2_PFC_hpc); cell2mat(cellfun(@(equis) length(unique(equis)),Cohfos2_PFC_hpc,'UniformOutput',false)) ])]; 
    
    TT.Properties.VariableNames=['Metric' cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)];
    
                writematrix(rat_no,strcat('ripple_spindle_HFO','.xls'),'Sheet',8,'Range',excel_range(2).(rat_no))
                writematrix('HPC Ripple - Spindles',strcat('ripple_spindle_HFO','.xls'), 'sheet',8, 'Range', 'B1')
                writetable(TT,strcat('ripple_spindle_HFO','.xls'),'Sheet',8,'Range',excel_range(1).(rat_no))    
                
%% HPC ripples-spindles-multiplets
    TT=table;
TT.Variables= [[{'PAR - Single'};{'PAR - doublet'};{'PAR - triplet'};{'PAR - quadruplet '}; {'PAR- Pentuplets'};{'PAR-Sextuplets'}; {' '}; {'PFC - Single'};{'PFC - doublet'};{'PFC - triplet'};{'PFC - quadruplet '}; {'PFC-pentplets'};{'PFC-Sextuplets'}], ...
                 num2cell([ [Out_PAR{1}(:,2) Out_PAR{2}(:,2) Out_PAR{3}(:,2) Out_PAR{4}(:,2); [NaN NaN NaN NaN] ;Out_PFC{1}(:,2) Out_PFC{2}(:,2) Out_PFC{3}(:,2) Out_PFC{4}(:,2)] ])];
    TT.Properties.VariableNames=['Metric' cellfun(@(equis) strrep(equis,'_','-'),g,'UniformOutput',false)];
    
                writematrix(rat_no,strcat('ripple_spindle_HFO','.xls'),'Sheet',9,'Range',excel_range(4).(rat_no))
                 writematrix('HPC Ripples-spindles multiplets',strcat('ripple_spindle_HFO','.xls'), 'sheet',9, 'Range', 'B1')
                writetable(TT,strcat('ripple_spindle_HFO','.xls'),'Sheet',9,'Range',excel_range(3).(rat_no))    

