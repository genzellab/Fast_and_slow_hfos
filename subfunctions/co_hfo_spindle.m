function [co_vec1,co_vec2]=co_hfo_spindle(a_s,a_m,a_e,N_s,N_m,N_e)%HPC,Cortex
    co_vec1=[];%HPC
    co_vec2=[];%Cortex
    parfor index_hfo=1:length(N_s)
     n=[N_s(index_hfo) N_e(index_hfo) ];   
     interval1 = fixed.Interval(n(1), n(2), '[]');
     for ind=1:length(a_s)
         
              interval2 = fixed.Interval(a_s(ind), a_e(ind), '[]');
              %Detect coocurrence if durations overlap
              if overlaps(interval1,interval2)
                co_vec1=[co_vec1 a_m(ind)];%HPC
                co_vec2=[co_vec2 N_m(index_hfo)];%Cortex

              end

     end
     
    end
end
