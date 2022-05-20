
function [sig]=getsignal(Sx,Ex,ti,V,k)
if ~isempty(Sx{k})
    for j=1:length(Sx{k})
    [~,ts]=min(abs(ti{k}-Sx{k}(j)));
    [~,tend]=min(abs(ti{k}-Ex{k}(j)));
    sig{j}=V{k}(ts:tend);
    end
else
    sig=[];
end
end







