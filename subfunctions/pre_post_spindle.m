function [Sx_pre,Mx_pre,Ex_pre,Sx_post,Mx_post,Ex_post]=pre_post_spindle(Sx,Mx,Ex)
 
dura=Ex-Sx; %duration

%Substract duration
Sx_pre=Sx-dura;
Mx_pre=Mx-dura;
Ex_pre=Ex-dura;

%Add duration
Sx_post=Sx+dura;
Mx_post=Mx+dura;
Ex_post=Ex+dura;

end