function [ Vg ] = bincube( v,ax )
%BINCUBE Summary of this function goes here
%   linear interpolation inside a cell (v)
%   Vg is the coordinates associated with the local index (ax)


s=ax(1);
t=ax(2);
u=ax(3);


x1=v(1,1);x2=v(1,2);x3=v(1,3);x4=v(1,4);x5=v(1,5);x6=v(1,6);x7=v(1,7);x8=v(1,8);
y1=v(2,1);y2=v(2,2);y3=v(2,3);y4=v(2,4);y5=v(2,5);y6=v(2,6);y7=v(2,7);y8=v(2,8);
z1=v(3,1);z2=v(3,2);z3=v(3,3);z4=v(3,4);z5=v(3,5);z6=v(3,6);z7=v(3,7);z8=v(3,8);

Vx=x1*(1-s)*(1-t)*(1-u)+x2*s*(1-t)*(1-u)+x3*s*(1-t)*u+x4*(1-s)*(1-t)*u+...
   x5*(1-s)*t*(1-u)+x6*s*t*(1-u)+x7*s*t*u+x8*(1-s)*t*u;
Vy=y1*(1-s)*(1-t)*(1-u)+y2*s*(1-t)*(1-u)+y3*s*(1-t)*u+y4*(1-s)*(1-t)*u+...
   y5*(1-s)*t*(1-u)+y6*s*t*(1-u)+y7*s*t*u+y8*(1-s)*t*u;
Vz=z1*(1-s)*(1-t)*(1-u)+z2*s*(1-t)*(1-u)+z3*s*(1-t)*u+z4*(1-s)*(1-t)*u+...
   z5*(1-s)*t*(1-u)+z6*s*t*(1-u)+z7*s*t*u+z8*(1-s)*t*u;
Vg=[Vx;Vy;Vz];





end

