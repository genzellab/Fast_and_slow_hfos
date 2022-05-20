function [channels]=gui_table_channels(z,rats,label1,str_table)
f = figure(1);
dat = z;
cnames=num2str(rats.');
cnames=strcat('Rat',cnames);
rnames=label1;

c = uicontrol('Style','text','Position',[1 380 300 20]);
c.String = {['Fill the table with the correct ' str_table]};
c.FontSize=10;
c.FontAngle='italic';
%%

t = uitable('Parent',f,'Data',dat,'ColumnName',cnames,... 
            'RowName',rnames,'Position',[20 280 450 100]);
set(t,'ColumnEditable',true(1,length(rats)))        


h = uicontrol('Position',[220 20 100 40],'String','Confirm',...
              'Callback','uiresume(gcbf)');
h.FontSize=10;
uiwait(gcf); 
aver= get(t,'Data');
aver(aver==0) = NaN;
close(f);


for i=1:size(cnames,1)
    cenames=cnames(i,:);
    cenames(isspace(cenames)) = [];
    Cnames{i,:}=cenames;
    
    channels(1).(Cnames{i})= aver(:,i).';
end
clear cenames

end