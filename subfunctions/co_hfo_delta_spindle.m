
function [co_vec1,co_vec2]=co_hfo_delta_spindle(a,N)%delta,spindle
    co_vec1=[];%delta
    co_vec2=[];%spindle
    for index_hfo=1:length(N);
     n=N(index_hfo);   

    [val,idx]=min(abs(a-n));
    minVal=a(idx);

    %Diference
    df=(minVal-n);
    %Coocur if within -0.5 to 1 sec difference
    if  df<=1 & df>=-0.50
        co_vec1=[co_vec1 minVal];
        co_vec2=[co_vec2 n];
    end

    end
end
