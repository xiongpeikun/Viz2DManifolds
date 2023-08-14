function [ stu ] = xyztostu( Nxyz,x,y,z,ijk,xyz )
%XYZTOSTU Summary of this function goes here
%   Detailed explanation goes here

r=Nxyz(4);

i=ijk(1);j=ijk(2);k=ijk(3);
x1=xyz(1,i,j,k);y1=xyz(2,i,j,k);z1=xyz(3,i,j,k);
s=(x-x1)/r;
t=(y-y1)/r;
u=(z-z1)/r;
stu=[s,t,u];

end

