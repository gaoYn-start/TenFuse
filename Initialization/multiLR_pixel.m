function [p]=multiLR_pixel(m,h,lbd1)
%%% m is column vector(4bands);
%%% h is a scale;


cons = h./(lbd1-m'*m);
p = cons*m;


 
 end
 
 
 
 
