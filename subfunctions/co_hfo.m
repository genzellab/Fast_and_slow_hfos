
function [co_vec1,co_vec2]=co_hfo(a,N)%HPC,Cortex
    co_vec1=[];%HPC
    co_vec2=[];%Cortex
    for index_hfo=1:length(N);
     n=N(index_hfo);   

    [val,idx]=min(abs(a-n));
    minVal=a(idx);

    %Diference
    df=abs(minVal-n);
    %Coocur if closer to 50ms
    if df<=0.050
        co_vec1=[co_vec1 minVal];
        co_vec2=[co_vec2 n];
    end

    end
end
