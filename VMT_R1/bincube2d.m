function [ Vg ] = bincube2d( v,ax )
%BINCUBE2D Summary of this function goes here
%   Detailed explanation goes here


s=ax(1);
t=ax(2);
u=ax(3);

x1=v(1);x2=v(2);x3=v(3);x4=v(4);x5=v(5);x6=v(6);x7=v(7);x8=v(8);
Vg=x1*(1-s)*(1-t)*(1-u)+x2*s*(1-t)*(1-u)+x3*s*(1-t)*u+x4*(1-s)*(1-t)*u+...
   x5*(1-s)*t*(1-u)+x6*s*t*(1-u)+x7*s*t*u+x8*(1-s)*t*u;



end

