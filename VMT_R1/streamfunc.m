function [  dy ] = streamfunc( t,stu,ijk,btotal)
%streamfunc Summary of this function goes here
%   Detailed explanation goes here

cell= asigncube( btotal, ijk );
dy=bincube( cell, stu);
dy=dy';


end

