function [ D ] = getjac(btotal,i,j,k)
%GETJAC Summary of this function goes here
%   Detailed explanation goes here


D(1,1)=(btotal(1,i+1,j,k)-btotal(1,i-1,j,k))/2;
D(1,2)=(btotal(1,i,j+1,k)-btotal(1,i,j-1,k))/2;
D(1,3)=(btotal(1,i,j,k+1)-btotal(1,i,j,k-1))/2;
D(2,1)=(btotal(2,i+1,j,k)-btotal(2,i-1,j,k))/2;
D(2,2)=(btotal(2,i,j+1,k)-btotal(2,i,j-1,k))/2;
D(2,3)=(btotal(2,i,j,k+1)-btotal(2,i,j,k-1))/2;
D(3,1)=(btotal(3,i+1,j,k)-btotal(3,i-1,j,k))/2;
D(3,2)=(btotal(3,i,j+1,k)-btotal(3,i,j-1,k))/2;
D(3,3)=(btotal(3,i,j,k+1)-btotal(3,i,j,k-1))/2;



end

