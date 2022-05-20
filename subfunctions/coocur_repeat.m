function [na1,na2]=coocur_repeat(a1,a2)

if length(unique(a1))~=length(a2) %Means a1 contains repeated values
    [na1,f2]=unique(a1);
    na2=[];
    for j=1:length(na1)
        [~,f2]=min(abs(a2-na1(j)));
        na2=[na2 a2(f2)];
    end
end

if length(unique(a2))~=length(a1)%Means a2 contains repeated values
    [na2,f2]=unique(a2);
    na1=[];
    for j=1:length(na2)
        [~,f2]=min(abs(a1-na2(j)));
        na1=[na1 a1(f2)];
    end
end

if length(unique(a2))==length(unique(a1)) %No repeated values or same amount of repeated values
    if  length(unique(a1))==length(a2) && length(unique(a2))==length(a1)
    na1=a1;
    na2=a2;
    else
            [na1,~]=unique(a1);
            [na2,~]=unique(a2);

        
    end
end

end