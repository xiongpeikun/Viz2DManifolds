function [ next,stu,ijk ] = insidecell( Nxyz,point,rv,rijk,xyz )
%INSIDECELL Summary of this function goes here
%   

Nx=Nxyz(1);
Ny=Nxyz(2);
Nz=Nxyz(3);
r=Nxyz(4);


next=0;
ijk=rijk;
stu=zeros(1,3);
x1=rv(1,1);x7=rv(1,7);
y1=rv(2,1);y7=rv(2,7);
z1=rv(3,1);z7=rv(3,7);

if point(1)<x1
    n=1+floor((x1-point(1))/r);
    ijk(1)=rijk(1)-n;
    next=1;
end
if point(1)>x7
    n=1+floor((point(1)-x7)/r);
    ijk(1)=rijk(1)+n;
    next=1;
end
if point(2)<y1
    n=1+floor((y1-point(2))/r);
    ijk(2)=rijk(2)-n;
    next=1;
end
if point(2)>y7
    n=1+floor((point(2)-y7)/r);
    ijk(2)=rijk(2)+n;
    next=1;
end
if point(3)<z1
    n=1+floor((z1-point(3))/r);
    ijk(3)=rijk(3)-n;
    next=1;
end
if point(3)>z7
    n=1+floor((point(3)-z7)/r);
    ijk(3)=rijk(3)+n;
    next=1;
end

if next==1
    if ijk(1)>Nx || ijk(1)<1 || ijk(2)>Ny||ijk(2)<1 || ijk(3)>Nz||ijk(3)<1
        return;
    end
end

i=ijk(1);j=ijk(2);k=ijk(3);
x1=xyz(1,i,j,k);y1=xyz(2,i,j,k);z1=xyz(3,i,j,k);
s=(point(1)-x1)/r;
t=(point(2)-y1)/r;
u=(point(3)-z1)/r;
stu=[s,t,u];

end

