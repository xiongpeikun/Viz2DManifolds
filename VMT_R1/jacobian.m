function [ j ] = jacobian( stu ,Ds )
%JACOBIAN Summary of this function goes here
%   Jacobian matrix

s=stu(1);t=stu(2);u=stu(3);

p1=(1-s)*(1-t)*(1-u);p2=s*(1-t)*(1-u);p3=s*(1-t)*u;p4=(1-s)*(1-t)*u;
p5=(1-s)*t*(1-u);p6=s*t*(1-u);p7=s*t*u;p8=(1-s)*t*u;

j=p1.*Ds(:,:,1)+p2.*Ds(:,:,2)+p3.*Ds(:,:,3)+p4.*Ds(:,:,5)+...
    p5.*Ds(:,:,5)+p6.*Ds(:,:,6)+p7.*Ds(:,:,7)+p8.*Ds(:,:,8);

end

