function [g1,g1_f,G,G_f,FB,FB1,n]=getval_granger(P,Q,labelconditions3,label1,s,w,fn)

%Find minimum number of events
    n=min([length(P.(labelconditions3{1}).(label1{w}){s}) length(P.(labelconditions3{2}).(label1{w}){s})...
        length(P.(labelconditions3{3}).(label1{w}){s}) length(P.(labelconditions3{4}).(label1{w}){s})]);
    


type_hfo{1}='cohfos';
type_hfo{2}='single';


for condition=1:length(labelconditions3)

    %Order ripples
    p=P.(labelconditions3{condition}).(label1{w}){s}; 
    q=Q.(labelconditions3{condition}).(label1{w}){s}; 
    % R=(cellfun(@(equis1) max(abs(hilbert(equis1(1,:)))),q));
    %R=(cellfun(@(equis1) max(abs(hilbert(equis1(1,121-50:121+50)))),q));
    R=(cellfun(@(equis1) max(abs(hilbert(equis1(1,151-25:151+25)))),q));
    [~,r]=sort(R,'descend');
    p=p(r);
    q=q(r);
    p=p(1:n);
    q=q(1:n);


    %Max 1000 ripples.
    if length(q)>1000
        q=q(1:1000);
        p=p(1:1000);
    end

    if w==1 %HPC-centered ripples.
         p=cellfun(@(equis1) flip(equis1),p,'UniformOutput',false);
         q=cellfun(@(equis1) flip(equis1),q,'UniformOutput',false);
    end
    
    ro=1200;
    %Compute non-parametric, parametric and conditional spectral GC.
    [gran,gran1,grangercon]=getgranger(p,create_timecell(ro,length(p),fn),'Wideband',ro,10,[0:1:300],fn);
    
    %Extract granger spectrums
    G{condition}=gran.grangerspctrm;%Non-parametric (Pairwise)
    g1{condition}=gran1.grangerspctrm;%Parametric
    g2{condition}=grangercon.grangerspctrm;%Non-parametric (Conditional)
    
 %Compute mean granger for different frequency bands.   
[FB{condition}]=gc_freqbands(gran,0,'granger');%Non-parametric (Pairwise)
[FB1{condition}]=gc_freqbands(gran1,0,'granger');%Parametric

    
end
%Frequencies _f
G_f=gran.freq;
g1_f=gran1.freq;
                  
end
