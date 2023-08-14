function [ x,y,yy,ijk,originflg,boundaryflg ] = rk4( Nx,Ny,Nz,q,h,ti,ijk,xyz,btotal )
%RK4 Summary of this function goes here
%   4th order Runge-Kutta method to integrate flow orbits
%   

originflg=0;
boundaryflg=0;
back=0;
cord= asigncube( xyz, ijk );
a1=ti(1);b1=ti(2);
if b1<0
   b1=-b1;
   back=1;
end
n=floor((b1-a1)/h); %step
x=zeros(1,n+1);
y=zeros(n+1,3);
yy=zeros(n+1,3);
x(1)=a1;    %init time
y(1,:)=q; %init value
yy(1,:)=bincube( cord, q );

for ii=1:n
    
    x(ii+1)=x(ii)+h;
    
    k1=streamfunc(x(ii),y(ii,:),ijk,btotal );
    
    k2=streamfunc(x(ii)+h/2,y(ii,:)+h*k1/2,ijk,btotal);
    
    k3=streamfunc(x(ii)+h/2,y(ii,:)+h*k2/2,ijk,btotal);
    
    k4=streamfunc(x(ii)+h,y(ii,:)+h*k3,ijk,btotal);
    if back==1
        k1=-k1;
        k2=-k2;
        k3=-k3;
        k4=-k4;
    end
    
    y(ii+1,:)=y(ii,:)+h*(k1+2*k2+2*k3+k4)/6;
    
    [ next,y(ii+1,:),ijk ] = whichcell( y(ii+1,:),ijk );
    
    if next==0
        yy(ii+1,:)=bincube( cord, y(ii+1,:) );
        if yy(ii+1,1)==yy(ii,1)&&yy(ii+1,2)==yy(ii,2)&&yy(ii+1,3)==yy(ii,3)
            x=x(1:ii);
            y=y(1:ii,:);
            yy=yy(1:ii,:);
            originflg=1;
            return;
        end
    end
    
    if next==1
        if ijk(1)>Nx || ijk(1)<1 || ijk(2)>Ny||ijk(2)<1 || ijk(3)>Nz||ijk(3)<1
            x=x(1:ii);
            y=y(1:ii,:);
            yy=yy(1:ii,:);
            boundaryflg=1;
            return;
        end
        cord= asigncube( xyz, ijk );
        yy(ii+1,:)=bincube( cord, y(ii+1,:) );
    end

    
end

end

