
function [C]=create_timecell(ro,leng,fn)
%create_timecell(ro,leng)
%iNPUTS:
%ro:1200
%leng:length(p)
    %fn=1000;
    vec=-ro/fn:1/fn:ro/fn;
    C    = cell(1, leng);
    C(:) = {vec};
end


