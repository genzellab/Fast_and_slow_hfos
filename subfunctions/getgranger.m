
function [granger,granger1,granger_cond,granger_cond_multi]=getgranger(q,timecell,label,ro,ord,freqrange,fn)
%Computes multiple types of granger causality. 
%Mainly parametric and Non parametric.

data1.trial=q;
data1.time= timecell; 
data1.fsample=fn;
data1.label=cell(3,1);

data1.label{1}='PAR';
data1.label{2}='PFC';
data1.label{3}='HPC';

%Parametric model
 [granger1]=createauto(data1,ord,'yes');  
 
%Non parametric
[granger]=createauto_np(data1,freqrange,[]);

%Non parametric Conditional
[granger_cond]=createauto_np(data1,freqrange,'yes');

%Parametric with multivariate setting set on (testing purposes)
[granger_cond_multi]=createauto_cond_multivariate(data1,ord);

end