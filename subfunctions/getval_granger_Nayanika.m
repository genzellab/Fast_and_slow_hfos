function [n, varargout]=getval_granger_Nayanika(P,Q,labelconditions3,label1,s,w,fn,tf)
% g1,g1_f,G,G_f,FB,FB1,n, granger_tf
%Find minimum number of events
if tf==0
    n=min([length(P.all.(labelconditions3{1}).(label1{w})) length(P.all.(labelconditions3{2}).(label1{w}))...
        length(P.all.(labelconditions3{3}).(label1{w})) length(P.all.(labelconditions3{4}).(label1{w}))]);
else
   n=length( [P.all.(labelconditions3{1}).(label1{w}) ...
       P.all.(labelconditions3{2}).(label1{w})...
       P.all.(labelconditions3{3}).(label1{w}) ...
       P.all.(labelconditions3{4}).(label1{w}) ])
end

% 
% type_hfo{1}='cohfos';
% type_hfo{2}='single';


if tf==1
    
   p=[P.all.(labelconditions3{1}).(label1{w}) ...
            P.all.(labelconditions3{2}).(label1{w})...
            P.all.(labelconditions3{3}).(label1{w}) ...
            P.all.(labelconditions3{4}).(label1{w}) ]; 
   
   q=[Q.all.(labelconditions3{1}).(label1{w}) ...
            Q.all.(labelconditions3{2}).(label1{w})...
            Q.all.(labelconditions3{3}).(label1{w}) ...
            Q.all.(labelconditions3{4}).(label1{w}) ];  
   
      
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
   
   % selecting random m samples to compute time frequency gc such that m=50% of n.
   % iterate 15 times
   iter=30
   m=0.2*min(n,1000)
   grangerspctrm_concat = zeros(3,3,300,201,iter);
   for i=1:iter
       i
       randorder = randperm(length(q));
       temp_q = q(randorder);
       temp_p = p(randorder);
       q_= temp_q(1:m);
       p_=temp_p(1:m);
       
       if w==1 %HPC-centered ripples.
           p_=cellfun(@(equis1) flip(equis1),p_,'UniformOutput',false);
           q_=cellfun(@(equis1) flip(equis1),q_,'UniformOutput',false);
       end
       
       % Compute time frequency gc
       ro=1200;
       [granger_tf]=getgranger_tf(p_,create_timecell(ro,length(p_),fn),'Wideband',ro,10,[0:1:300],fn);
       
       grangerspctrm_concat(:,:,:,:,i)=granger_tf.grangerspctrm;
       
   end
   
   granger_tf.grangerspctrm = grangerspctrm_concat;
    
else
    % tf=0, Compute non-parametric, parametric and conditional spectral GC.
    for condition=1:length(labelconditions3)

        %Order ripples
        p=P.all.(labelconditions3{condition}).(label1{w}); 
        q=Q.all.(labelconditions3{condition}).(label1{w}); 
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
        [gran,gran1,grangercon, ~]=getgranger(p,create_timecell(ro,length(p),fn),'Wideband',ro,10,[0:1:300],fn);

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

if tf==1
       varargout{1}=granger_tf;
else
    
    varargout{1}=g1;
    varargout{2}=g1_f;
    varargout{3}=G;
    varargout{4}=G_f;
    varargout{5}=FB;
    varargout{6}=FB1;
end

end
