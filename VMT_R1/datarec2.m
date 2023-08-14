function [  ] = datarec2( arclength,adcr,cr,cbnum,ijkindex,btotal,xyz,p,moment,current,rho )
%DATAREC Summary of this function goes here
%   Detailed explanation goes here
Nxyz=[360,240,240,0.25];
% adijkindex=newijkindex;
% adstu=zero(numcb,3);
% newstu=adstu;
cr=cr';
lie=zeros(cbnum,14);
datainer=[adcr,lie];
dataouter=[cr,lie];
rannum=cbnum*2;
rangeddata=zeros(rannum,17); %"x","y","z","btotal","bx","by","bz","p","mvx"
                            %"mvy","mvz","mvtotal","jx","jy","jz","jtotal","n"

num=0;
for jj=1:cbnum
    rijk=ijkindex(jj,:);
    newcube=asigncube( xyz,rijk );
    [ ~,stu,nijk ] = insidecell( Nxyz,adcr(jj,:),newcube,rijk,xyz );
%     adijkindex(jj,:)=nijk;
%     adstu(jj,:)=stu;
    nbcube=asigncube(btotal,nijk);%btotalcube
    nb=bincube( nbcube,stu );
    npcube=asigncube2d(p,nijk);%pcube
    npb=bincube2d(npcube,stu);
    datainer(jj,4)=norm(nb); %nb(1)^2+nb(2)^2+nb(3)^2;  %btotal of this point
    datainer(jj,5)=nb(1);   %u
    datainer(jj,6)=nb(2);   %v
    datainer(jj,7)=nb(3);   %w
    datainer(jj,8)=npb;      %p
    
    nvcube=asigncube(moment,nijk);%mvcube
    nv=bincube(nvcube,stu);
    datainer(jj,9)=nv(1);      %mvx
    datainer(jj,10)=nv(2);      %mvy
    datainer(jj,11)=nv(3);      %mvz
    datainer(jj,12)=norm(nv);%nv(1)^2+nv(2)^2+nv(3)^2;      %mvtotal
    
    njcube=asigncube(current,nijk); %jcube
    nj=bincube( njcube,stu );
     datainer(jj,13)=nj(1);      %jx
    datainer(jj,14)=nj(2);      %jy
    datainer(jj,15)=nj(3);      %jz
    datainer(jj,16)=norm(nj);%nj(1)^2+nj(2)^2+nj(3)^2;      %jtotal
    
    nncube=asigncube2d(rho,nijk);%ncube
    nnb=bincube2d(nncube,stu);
    datainer(jj,17)=nnb;
    
    nstu = xyztostu( Nxyz,cr(jj,1),cr(jj,2),cr(jj,3),rijk,xyz );
%     newstu(jj,:)=nstu;
    bcube=asigncube(btotal,rijk);
    cbb=bincube(bcube,nstu);
    pcube=asigncube2d(p,rijk);
    pb=bincube2d(pcube,nstu);
    dataouter(jj,4)=norm(cbb);%cbb(1)^2+cbb(2)^2+cbb(3)^2;
    dataouter(jj,5)=cbb(1);
    dataouter(jj,6)=cbb(2);
    dataouter(jj,7)=cbb(3);
    dataouter(jj,8)=pb;
    
    vcube=asigncube(moment,rijk);%mvcube
    pv=bincube(vcube,nstu);
    dataouter(jj,9)=pv(1);
    dataouter(jj,10)=pv(2);
    dataouter(jj,11)=pv(3);
    dataouter(jj,12)=norm(pv);%pv(1)^2+pv(2)^2+pv(3)^2;
    
    jcube=asigncube(current,rijk); %jcube
    j=bincube( jcube,nstu );
    dataouter(jj,13)=j(1);     
    dataouter(jj,14)=j(2);     
    dataouter(jj,15)=j(3);      
    dataouter(jj,16)=norm(j);%j(1)^2+j(2)^2+j(3)^2;      
    
    ncube=asigncube2d(rho,rijk);%ncube
    n=bincube2d(ncube,nstu);
    dataouter(jj,17)=n;
    
    num=num+1;
    rangeddata(num,:)=datainer(jj,:);
    num=num+1;
    rangeddata(num,:)=dataouter(jj,:);
    
    
end

ring=[rangeddata;rangeddata(1:2,:)];


% header = ['x','y','z','btotal'];
hou='.dat';
name=num2str(arclength);
path='.\plotfile2\';
namef=[path,name,hou];
fid = fopen(namef,'w');

fprintf(fid,'variables="x","y","z","btotal","u","v","w","p","mvx","mvy","mvz","mvtotal","jx","jy","jz","jtotal","n"\n');
numk=cbnum+1;
fprintf(fid,'zone t="Frame %s" i=2,j=%d,k=1\n',name,numk);

for m = 1:size(ring,1)
  fprintf(fid,'%f %f %f %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n',ring(m,1),ring(m,2),ring(m,3),ring(m,4),ring(m,5),ring(m,6),ring(m,7),ring(m,8),ring(m,9),ring(m,10),ring(m,11),ring(m,12),ring(m,13),ring(m,14),ring(m,15),ring(m,16),ring(m,17));
end
fclose(fid);



end

