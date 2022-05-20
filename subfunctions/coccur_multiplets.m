function [out_par]=coccur_multiplets(cohfos2_hpc)

    rep_coccur=diff(sort([cohfos2_hpc{:}]));

    rep_coccur(rep_coccur~=0)=1;

    rep_coccur=not(rep_coccur);

     rep_coccur=ConsecutiveOnes(rep_coccur);
    rep_coccur=rep_coccur(rep_coccur~=0)+1 ;%Ignore single values and add one to label as number of coocur

    a = 1:6;
    hi=histc(rep_coccur(:),a);
    if size(hi,2)> size(hi,1)
        hi=hi.';
    end
    out_par = [a.',hi];
    out_par(1,2)=length([cohfos2_hpc{:}])-sum(out_par(:,1).*out_par(:,2)); %Add the extra coocur lost when using diff()

end