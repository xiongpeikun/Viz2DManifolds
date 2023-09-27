function [points, brs,oflg] = intersectfr( Nxyz,qstu,q,r,fn,ti,h,ijk,delta,arctheta,dirbr,xyz,btotal )
%INTERSECTFR Summary of this function goes here
%   calculate all intersection of a flow form q
%   q is the init point of flow



Nx=Nxyz(1);
Ny=Nxyz(2);
Nz=Nxyz(3);
mex=Nxyz(5);

e=0.01;
deltam=delta*(1+e);
origin=[0,0,0];

eps=3e-14;
oq=q;
oh=h;
oti=ti;


grow=0.1+3*norm(q-r)+delta*5; % computation distance of flow
fnt=fn';
mo=norm(fn);
points=0;
dis=0;
% dis2=0;
brs=[];
oflg=0;
stopf=0;

% intersections of a flow
step=0;
while dis<=grow&&stopf==0
    if step==60000
        break;
    end
    if (ijk(1)>Nx||ijk(1)<1) || (ijk(2)>Ny||ijk(2)<1) || (ijk(3)>Nz||ijk(3)<1)
        break;
    end
    
    if norm(origin-q)>3
        [ti,h,~] = timestp( Nxyz,q,qstu,ti,h,ijk,delta,xyz,btotal );
    end
    ijkb=ijk;
    if mex==1
        [ ~,a,aa,ijk,oflg,bflg] = rk4_mex( Nx,Ny,Nz,qstu,h,ti,ijk,xyz,btotal );
    else
        [ ~,a,aa,ijk,oflg,bflg] = rk4( Nx,Ny,Nz,qstu,h,ti,ijk,xyz,btotal );
    end
    dislen=norm(aa(1,:)-aa(end,:));
    if dislen>delta*3
        n=dislen/delta;
        tspan=ti(2)/n;
        h=h/n;
        ti=[0,tspan];
        if mex==1
            [ ~,a,aa,ijk,oflg,bflg] = rk4_mex( Nx,Ny,Nz,qstu,h,ti,ijkb,xyz,btotal );
        else
            [ ~,a,aa,ijk,oflg,bflg] = rk4( Nx,Ny,Nz,qstu,h,ti,ijkb,xyz,btotal );
        end
    end
    
%     plot3(aa(:,1),aa(:,2),aa(:,3));
%     hold on;

    h=oh;
    ti=oti;
    
    qstu=a(end,:);
    q=aa(end,:);
    dis=norm(oq-q);
    
    id=length(aa(:,1));
    ptp=zeros(1,id);
    
    for i=1:id
        ptp(i)=(r-aa(i,:))*fnt/mo; %the destances between points and Fr plane
    end
    
    if id==length(find(ptp<-eps))||id==length(find(ptp>eps))% for no intersection
        if norm(aa(1,:)-q)<1e-15
            break;
        end
        step=step+1;
        if oflg==1||bflg==1
            break;
        end
        continue;
    end
    
    [numpts,interpts] = shootm( ptp,aa,id,r,fn ); %looking for intersection
    
    
    for i=1:numpts
         tmpdir=interpts(i,:)-r;
         disarc=norm(tmpdir);
        if ((dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir))<arctheta)
            if disarc>deltam
                stopf=1;
            end
        end

    end

    points=points+numpts;
    brs=cat(1,brs,interpts);
    
    if norm(q-origin)<=3
        break;
    end
    if norm(aa(1,:)-q)<1e-15
        break;
    end
    step=step+1;
    if oflg==1||bflg==1
        break;
    end
    
end

if points>1
    pnum=0;
    inbr=zeros(points,3);
    brs=[brs;brs(1,:)];
    for i=2:points
        if norm(brs(i-1,:)-brs(i,:))<eps*10
            brs(i,:)=brs(i-1,:);
            if i==points
                pnum=pnum+1;
                inbr(pnum,:)=brs(i,:);
            end
            continue;
        end
        pnum=pnum+1;
        inbr(pnum,:)=brs(i-1,:);
        if i==points
            pnum=pnum+1;
            inbr(pnum,:)=brs(i,:);
        end
    end
    inbr=inbr(1:pnum,:);
    points=pnum;
    brs=inbr;
end




end

