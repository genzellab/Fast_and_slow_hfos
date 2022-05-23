

function granger_plot(g,g_f,labelconditions,freqrange)
%Plots granger values across frequencies

allscreen()
myColorMap=StandardColors; 
F= [1 2; 1 3; 2 3] ;

%Labels
lab=cell(6,1);

lab{2}='PFC -> PAR';
lab{1}='PAR -> PFC';

lab{4}='HPC -> PAR';
lab{3}='PAR -> HPC';

lab{6}='HPC -> PFC';
lab{5}='PFC -> HPC';
% 
 for j=1:3

     f=F(j,:);

     mmax1=max([max(squeeze(g{1}(f(1),f(2),:))) max(squeeze(g{2}(f(1),f(2),:))) ...
         max(squeeze(g{3}(f(1),f(2),:))) max(squeeze(g{4}(f(1),f(2),:)))]);

     mmax2=max([max(squeeze(g{1}(f(2),f(1),:))) max(squeeze(g{2}(f(2),f(1),:))) ...
         max(squeeze(g{3}(f(2),f(1),:))) max(squeeze(g{4}(f(2),f(1),:)))]);

     mmax=max([mmax1 mmax2]);


     subplot(3,2,2*j-1)
     plot(g_f, squeeze(g{1}(f(1),f(2),:)),'LineWidth',2,'Color',myColorMap(1,:))
     hold on
     plot(g_f, squeeze(g{2}(f(1),f(2),:)),'LineWidth',2,'Color',myColorMap(2,:))
     plot(g_f, squeeze(g{3}(f(1),f(2),:)),'LineWidth',2,'Color',myColorMap(3,:))
     plot(g_f, squeeze(g{4}(f(1),f(2),:)),'LineWidth',2,'Color',myColorMap(4,:))

     xlim(freqrange)

     xlabel('Frequency (Hz)')
     ylabel('G-causality')
     title(lab{2*j-1})
     subplot(3,2,2*j)
     plot(g_f, squeeze(g{1}(f(2),f(1),:)),'LineWidth',2,'Color',myColorMap(1,:))
     hold on
     plot(g_f, squeeze(g{2}(f(2),f(1),:)),'LineWidth',2,'Color',myColorMap(2,:))
     plot(g_f, squeeze(g{3}(f(2),f(1),:)),'LineWidth',2,'Color',myColorMap(3,:))
     plot(g_f, squeeze(g{4}(f(2),f(1),:)),'LineWidth',2,'Color',myColorMap(4,:))


     xlim(freqrange)
     xlabel('Frequency (Hz)')
     ylabel('G-causality')
    if j==1
     legend(labelconditions,'Location','best') %Might have to change to default. 
    end
    title(lab{2*j})

 end
 
end