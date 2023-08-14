function [ ti,h,con] = timestp( Nxyz,r,rstu,ti,h,rijk,delta,xyz,btotal)
%TIMESTP Summary of this function goes here
%   refine time step
% 

Nx=Nxyz(1);
Ny=Nxyz(2);
Nz=Nxyz(3);

con=0;
back=0;

tspan=ti(1,end);

bu=streamfunc( 1,rstu,rijk,btotal);
if abs(bu(1))>1e-8||abs(bu(2))>1e-8||abs(bu(3))>1e-8
    h=h/2;
    tspan=tspan/2;
    ti=[0,tspan];
end
if tspan<0
    tspan=-tspan;
    back=1;
end
hh=h;
tti=ti;
[ ~,~,at,~,oflg] = rk4_mex( Nx,Ny,Nz,rstu,h,ti,rijk,xyz,btotal );
if oflg==1&&length(at(:,1))>1&&length(at(:,1))<41
    len=length(at(:,1))-1;
    tspan=tspan/40*len;
    
elseif length(at(:,1))==1
    con=1;
    return;
else
    flowlen=norm(at(end,:)-r);
    
    % plot3(at(:,1),at(:,2),at(:,3));
    % hold on;
    
    arclen=delta;
    n=flowlen/arclen;
    tspan=tspan/n;
    h=tspan/40;
end
if back==1
    tspan=-tspan;
end
ti=[0,tspan];
% testti=[0,h*2];
[ ~,~,at,~,oflg2] = rk4_mex( Nx,Ny,Nz,rstu,h,ti,rijk,xyz,btotal );
if length(at(:,1))==1&&oflg==oflg2&&oflg==1
    con=1;
    return;
end

if length(at(:,1))==1
    h=hh;
    ti=tti;
end
[ ~,~,at,~,~] = rk4_mex( Nx,Ny,Nz,rstu,h,ti,rijk,xyz,btotal );
dflow=norm(at(1,:)-at(2,:));

dmax=5e-5;
dmin=1e-5;

if dflow>dmax
    
    bmax=dflow/dmax;
    h=h/bmax;
%     if back==1
%         testti=[0,-h];
%     else
%         testti=[0,h];
%     end
%     h=h/2;
%     [ ~,~,at,~] = rk4_mex( Nx,Ny,Nz,rstu,h,testti,rijk,xyz,btotal );
%     if length(at(:,1))==1
%         break;
%     end
%     dflow=norm(at(1,:)-at(2,:));
    
end



while dflow<dmin&&dflow>0
    
    h=h*2;
    if back==1
        testti=[0,-2*h];
    else
        testti=[0,2*h];
    end
    [ ~,~,at,~] = rk4_mex(Nx,Ny,Nz,rstu,h,testti,rijk,xyz,btotal );
    if length(at(:,1))==1
        h=h/2;
        break;
    end
    dflow=norm(at(1,:)-at(2,:));
    
end



end

