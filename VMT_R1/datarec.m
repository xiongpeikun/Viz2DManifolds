function [  ] = datarec( Nxyz,arclength,newadcr,cbr,cbnum,newijkindex,btotal,xyz )
%DATAREC Summary of this function goes here
%   Generating binary file of Geodesic bands for TecPlot.

cbr=cbr';
lie=zeros(cbnum,4);
datainer=[newadcr,lie];
dataouter=[cbr,lie];
rannum=cbnum*2;
rangeddata=zeros(rannum,7);

num=0;
for jj=1:cbnum
    rijk=newijkindex(jj,:);
    newcube=asigncube( xyz,rijk );
    [ ~,stu,nijk ] = insidecell( Nxyz,newadcr(jj,:),newcube,rijk,xyz );
    nbcube=asigncube(btotal,nijk);
    nb=bincube( nbcube,stu );
    datainer(jj,4)=norm(nb);  %btotal of this point
    datainer(jj,5)=nb(1);   %u
    datainer(jj,6)=nb(2);   %v
    datainer(jj,7)=nb(3);   %w
    
    
    nstu = xyztostu( Nxyz,cbr(jj,1),cbr(jj,2),cbr(jj,3),rijk,xyz );
    bcube=asigncube(btotal,rijk);
    cbb=bincube(bcube,nstu);
    dataouter(jj,4)=norm(cbb);
    dataouter(jj,5)=cbb(1);
    dataouter(jj,6)=cbb(2);
    dataouter(jj,7)=cbb(3);
    
    num=num+1;
    rangeddata(num,:)=datainer(jj,:);
    num=num+1;
    rangeddata(num,:)=dataouter(jj,:);
    
    
end

ring=[rangeddata;rangeddata(1:2,:)];


hou='.dat';
name=num2str(arclength);
path='.\plotfile\';
namef=[path,name,hou];
fid = fopen(namef,'w');

fprintf(fid,'variables="x","y","z","btotal","u","v","w"\n');
numk=cbnum+1;
fprintf(fid,'zone t="Frame %s" i=2,j=%d,k=1\n',name,numk);

for m = 1:size(ring,1)
  fprintf(fid,'%f %f %f %e %e %e %e\n',ring(m,1),ring(m,2),ring(m,3),ring(m,4),ring(m,5),ring(m,6),ring(m,7));
end
fclose(fid);



end

