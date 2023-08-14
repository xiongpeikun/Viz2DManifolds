%Using for creating Mex file.
clc;clear;
load arclen;
name=num2str(arclength);
path='.\bands\';
hou='.mat';
nameb=[path,name,hou];
load('datas/271xyz');load('datas/271btotal');load (nameb);
numcr=length(cr);
tspan=-50000000;                 %time step of runge_kutta
ti=[0,tspan];               %init_time span
h=1000000;             %init time step
Nx=240;
Ny=160;
Nz=160;
r=0.25;
Nxyz=[Nx,Ny,Nz,r];
% delta=0.125;
crstu=zeros(numcr,3);
for ii=1:numcr
crstu(ii,:) = xyztostu( Nxyz,cr(1,ii),cr(2,ii),cr(3,ii),ijkindex(ii,:),xyz );
end
nijkindex=ijkindex;

% for ii=1:numcr%121:180
    ii=1;
    ijk=nijkindex(ii,:);% need to save every r 's ijk
    qstu=crstu(ii,:);
%     ti=[0,tspan];
%     h=1250000;
%     [ti,h,~] = timestp( Nxyz,cr(:,ii)',qstu,ti,h,ijk,delta,xyz,btotal );
    
%     for spi=1:100
%         if (ijk(1)>Nx||ijk(1)<1) || (ijk(2)>Ny||ijk(2)<1) || (ijk(3)>Nz||ijk(3)<1)
%             break;
%         end
%         [ti,h,~] = timestp( Nxyz,cr(:,ii)',qstu,ti,h,ijk,delta,xyz,btotal );
        [ ~,a,aa ,ijk] = rk4( Nx,Ny,Nz,qstu,h,ti,ijk,xyz,btotal );
%         plot3(aa(:,1),aa(:,2),aa(:,3));
%         hold on;
%         qstu=a(end,:);
% %         [ti,h,~] = timestp( Nxyz,a,qstu,ti,h,ijk,delta,xyz,btotal );
%     end

% end

