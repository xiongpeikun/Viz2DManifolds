% MIT License
% Copyright (c) 2023 Xiong Peikun

% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

% instruction for running code with Mex :
% Strongly recommend to compile the Mex file of rk4.m in Matlab Coder 
% Before compiling the Mex file, you must run initcircles.m to generate
% the necessary varibles.
% You can use flowline.m for Autodefine Input Types and checking issues.
% the name of Mex function should be rk4_mex (default name).
% Be aware you need to recompile and overwrite the mex file when you change 
% the size of cubic data (Nxyz).
%
% If you don't want to compile the Mex file. Please change "rk4_mex"
% function to "rk4" in files timestp.m and intersectfr.m

% Email: crazyxpk[at]gmail.com  please change [at] to @

clc;clear;

load('datas/271xyz');load('datas/271btotal');load('datas/nulldata271');
% xyz(3,x,y,z);  btotal(3,x,y,z);  nullpoints;nullijk;nullstus are (:,3)

Nx=length(xyz(1,:,1,1))-1;
Ny=length(xyz(1,1,:,1))-1;
Nz=length(xyz(1,1,1,:))-1;
rcube=norm(xyz(:,1,1,1)-xyz(:,1,1,2));  %cubic cell size
Nxyz=[Nx,Ny,Nz,rcube];

parnum=7; %number of workers for parfor function,
          %for non-parallel computation set parnum=1;          
deltao=0.5; %geodesic step size

cpi=1; %nullpoints index

%constructing initial circle
initcircles(cpi,Nxyz,deltao,nullijk,nullpoints,nullstus,xyz,btotal);

arclength=80; % the arc length tend to construct
orbitdflg=1; %boolean 1 for B(s) type CPs, 0 for A(s) type CPs

%constructing GLS manifold
leaf(arclength,Nxyz,orbitdflg,nullpoints,parnum,deltao,xyz,btotal);


