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
% Strongly recommend to compile the Mex file of rk4.m in Matlab Coder (using flowline.m)
% the name of Mex function should be "rk4_mex" (default name).
% If you don't want to compile the Mex file. Please set mex=0 below;
% 
% Email: crazyxpk[at]gmail.com  please change [at] to @

clc;clear;

load('datas/271xyz');load('datas/271btotal');load('datas/nulldata271');
% xyz(3,x,y,z);  btotal(3,x,y,z);  nullpoints;nullijk;nullstus are (:,3)

Nx=length(xyz(1,:,1,1))-1;
Ny=length(xyz(1,1,:,1))-1;
Nz=length(xyz(1,1,1,:))-1;
rcube=norm(xyz(:,1,1,1)-xyz(:,1,1,2));  %cubic cell size
mex=0; %set 1 to enable rk4_mex file, 0 to disable. 
Nxyz=[Nx,Ny,Nz,rcube,mex];

parnum=7; %number of workers for parfor function,
          %for non-parallel computation set parnum=1;
          
ir=0.5;  %initial circle radius          
deltao=0.5; %geodesic step size
trand=[deltao*0.5,deltao*1.5]; %transverse distance size of 
                                    % adjacent mesh points [min,max]
                                    
orbitdflg=0; %boolean 1 for B(s) type CPs, 0 for A(s) type CPs
cpi=3; %nullpoints index

tarclength=80; % the arc length tend to construct

%constructing initial circle
initcircles(cpi,Nxyz,ir,deltao,nullijk,nullpoints,nullstus,xyz,btotal);

%constructing GLS manifold
leaf(tarclength,Nxyz,orbitdflg,nullpoints,parnum,trand,deltao,xyz,btotal);



