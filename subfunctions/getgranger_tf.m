function [granger_tf]=getgranger_tf(q,timecell,label,ro,ord,freqrange,fn)
%Computes multiple types of granger causality. 
%Mainly parametric and Non parametric.

data1.trial=q;
data1.time= timecell; 
data1.fsample=fn;
data1.label=cell(3,1);

data1.label{1}='PAR';
data1.label{2}='PFC';
data1.label{3}='HPC';

% time frequency 
[granger_tf, freq_tf]=createauto_timefreq_Nayanika(data1,freqrange);
end

