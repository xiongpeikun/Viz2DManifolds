function [cube] = asigncube2d( data, a )
%ASIGNCUBE2D Summary of this function goes here
%   Detailed explanation goes here
cube=zeros(1,8);

i=a(1);j=a(2);k=a(3);

cube(:,1)=data(:,i,j,k);
cube(:,2)=data(:,i+1,j,k);
cube(:,3)=data(:,i+1,j,k+1);
cube(:,4)=data(:,i,j,k+1);
cube(:,5)=data(:,i,j+1,k);
cube(:,6)=data(:,i+1,j+1,k);
cube(:,7)=data(:,i+1,j+1,k+1);
cube(:,8)=data(:,i,j+1,k+1);

end

