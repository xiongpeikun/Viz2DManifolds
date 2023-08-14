
function [ ] = initcircles( spi,Nxyz,deltao,nullijk,nullpoints,nullstus,xyz,btotal )
%   Generate initial circle

% load('datas/184xyz_8');load('datas/184btotal_8');load('datas/nulldata184_8');

% deltao=0.5; %geodesic step
r=deltao;      %Radius

Ds=zeros(3,3,8);

i=nullijk(spi,1);j=nullijk(spi,2);k=nullijk(spi,3);
Ds(:,:,1)=getjac(btotal,i,j,k);
Ds(:,:,2)=getjac(btotal,i+1,j,k);
Ds(:,:,3)=getjac(btotal,i+1,j,k+1);
Ds(:,:,4)=getjac(btotal,i,j,k+1);
Ds(:,:,5)=getjac(btotal,i,j+1,k);
Ds(:,:,6)=getjac(btotal,i+1,j+1,k);
Ds(:,:,7)=getjac(btotal,i+1,j+1,k+1);
Ds(:,:,8)=getjac(btotal,i,j+1,k+1);
stu=nullstus(spi,:);
jtensor = jacobian( stu ,Ds );
% deter=det(jtensor);
[eigvec,eigval] = eig(jtensor);

if eigval(1)*eigval(5)>0
    evec1=eigvec(:,1);
    evec2=eigvec(:,2);
elseif eigval(1)*eigval(9)>0
    evec1=eigvec(:,1);
    evec2=eigvec(:,3);
else
    evec1=eigvec(:,2);
    evec2=eigvec(:,3);
end
origin=nullpoints(spi,:);

normal=cross(evec1,evec2);
% r=deltao;

ijk=nullijk(spi,:);
thiscord=asigncube( xyz, ijk );
% thiscell=asigncube( btotal, ijk );

u=[normal(2);-normal(1);0];
cx=normal(1)*normal(3);
cy=normal(2)*normal(3);
cz=-normal(1)*normal(1)-normal(2)*normal(2);
cxyz=[cx;cy;cz];

u=u/sqrt(normal(2)^2+normal(1)^2);
v=cxyz/sqrt(cx*cx+cy*cy+cz*cz);
n=21;
crpstu=zeros(3,n);
crp=zeros(3,n);
h=2*pi/(n-1);

for i=1:n
    t=(i-1)*h;
    crp(1,i)=origin(1)+r*(u(1)*cos(t)+v(1)*sin(t));
    crp(2,i)=origin(2)+r*(u(2)*cos(t)+v(2)*sin(t));
    crp(3,i)=origin(3)+r*(u(3)*cos(t)+v(3)*sin(t));
    [ ~,stu1,~ ] = insidecell( Nxyz,crp(:,i),thiscord,ijk,xyz );
    crpstu(:,i)=stu1;
end

plot3(crp(1,:),crp(2,:),crp(3,:),'color',[1 0 0]);% 1 0 0 ;0.5 0.5 1; 0 1 0
hold on;
xlabel('x');
ylabel('y');
zlabel('z');

cr=crp(:,1:end-1);
len=length(cr);
% crstu=crpstu(:,1:end-1)';
adcr=repmat(origin,len,1);
deltas=repmat(deltao,1,len);
ijkindex=repmat(ijk,len,1);
for jj=1:len
    rijk=ijkindex(jj,:);
    rcube=asigncube( xyz,rijk );
    [ ~,~,nijk ] = insidecell( Nxyz,cr(:,jj),rcube,rijk,xyz );
    ijkindex(jj,:)=nijk;
end

lost=zeros(1,len);
arclength=r;
cbnum=length(cr(1,:));

datarec(Nxyz,arclength,adcr,cr,cbnum,ijkindex,btotal,xyz);
arcthetas=repmat(0.2,1,cbnum);
polanums=zeros(3,cbnum);
name=num2str(arclength);
path='.\bands\';
hou='.mat';
nameb=[path,name,hou];
save(nameb,'cr','adcr','deltas','lost','ijkindex','arclength','arcthetas','polanums');
save('arclen','arclength');
    
end

