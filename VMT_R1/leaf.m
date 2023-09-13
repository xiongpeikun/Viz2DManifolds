function [ ] = leaf( tendarc,Nxyz,orbitdflg,nullpoints,parnum,trand,deltao,xyz,btotal )
%LEAF Summary of this function goes here
%   Main file of the GLS method

load arclen;
name=num2str(arclength);
path='.\bands\';
hou='.mat';
nameb=[path,name,hou];
load (nameb);
cr=cr;adcr=adcr;deltas=deltas;lost=lost;arclength=arclength;
ijkindex=ijkindex;arcthetas=arcthetas;polanums=polanums;

numcr=length(cr(1,:));
if orbitdflg==0
    tspan=-50000000;                 %init time step of runge_kutta
else
    tspan=50000000; 
end
ti=[0,tspan];               %init_time span
h=1250000;             %init time step
fprintf('arclength = %g\n ', arclength);
taomaxs=zeros(1,numcr);
% deltao=0.5;            %init arclength in every step
deltamin=deltao/4;         %minimum arclength
[ deltas,len,con,~] = adjdelta( cr,adcr,deltas,taomaxs,deltamin); 
if numcr==20
    deltas=repmat(deltamin,1,numcr);
end
deltalen=len;
% parnum=7;
hbar=waitbar(0,'initializing..');

polanums=zeros(3,numcr);

while arclength<tendarc
    
    tao=3;                 %number of points between rl and r 
    adelta=deltas;
    br=zeros(numcr,3);
    taomaxs=zeros(1,numcr);
    oijkindex=ijkindex;
    oarcthetas=arcthetas;
    opolanums=polanums;
    bvpflg=zeros(1,numcr);
    waitbar(0,hbar,'BVP is working...');
    ensti=num2str(numcr);
    
    for i=1:numcr
        
        sti=num2str(i);
        v=['BVP is working...',sti,'/',ensti];
        waitbar(i/numcr,hbar,v);
        
        if bvpflg(i)==1
            continue;
        end

        if con(i)==1 %if convergence detected on this leaf,then skip computation
            [ r,~,~ ] = selectr( i,numcr,cr,ijkindex ); 
            br(i,:)=r;
            bvpflg(i)=1;
            polanums(:,i)=zeros(3,1);
            continue;
        end
        
        parn=parnum+i-1;
        if parn>numcr
            parn=numcr;
        end
        if parnum>1
            parfor ii=i:parn %adding a new circle; parallel computarion available (parfor)
                if con(ii)==1 %if convergence detected on this leaf,then skip computation
                    [ r,~,~ ] = selectr( ii,numcr,cr,ijkindex );
                    br(ii,:)=r;
                    bvpflg(ii)=1;
                    polanums(:,ii)=zeros(3,1);
                    continue;
                end
                [accbr,taomax,newdelta,newarctheta,polanum] = leafbvp( Nxyz,h,ti,ii,numcr,cr,tao,adelta,deltao,adcr,deltamin,ijkindex,oarcthetas,nullpoints,opolanums,xyz,btotal ); %boundry value problem for computing br
                arcthetas(ii)=newarctheta;
                br(ii,:)=accbr;
                polanums(:,ii)=polanum;
                deltas(ii)=newdelta;
                taomaxs(ii)=taomax;
                bvpflg(ii)=1;
            end
        else
             for ii=i:parn %adding a new circle; parallel computarion available (parfor)
                if con(ii)==1 %if convergence detected on this leaf,then skip computation
                    [ r,~,~ ] = selectr( ii,numcr,cr,ijkindex );
                    br(ii,:)=r;
                    bvpflg(ii)=1;
                    polanums(:,ii)=zeros(3,1);
                    continue;
                end
                [accbr,taomax,newdelta,newarctheta,polanum] = leafbvp( Nxyz,h,ti,ii,numcr,cr,tao,adelta,deltao,adcr,deltamin,ijkindex,oarcthetas,nullpoints,opolanums,xyz,btotal ); %boundry value problem for computing br
                arcthetas(ii)=newarctheta;
                br(ii,:)=accbr;
                polanums(:,ii)=polanum;
                deltas(ii)=newdelta;
                taomaxs(ii)=taomax;
                bvpflg(ii)=1;
            end
        end

        if parn==numcr
            break;
        end
        
    end
        

    waitbar(1,hbar,'BVP is completed');
    
  
%     brr=[br;br(1,:)];
%     plot3(brr(:,1),brr(:,2),brr(:,3));
    if length(unique(deltas))>1
        fprintf('maxtao detected\n');
     end
    
    [newpolanums,cbr,newadcr,newdeltas,newtaomaxs,newlost,newijkindex,newarcthetas]= adjcbr_new( trand,Nxyz,br,parnum,cr,tao,ti,h,deltas,deltao,adcr,deltamin,taomaxs,lost,ijkindex,oijkindex,arcthetas,nullpoints,polanums,opolanums,con,xyz,btotal);% mesh adaptation
     
    plot3(cbr(:,1),cbr(:,2),cbr(:,3),'color',[1 0 0]); %plot the new circle
    drawnow;
    
    cbnum=length(cbr(:,1))-1;
    
    if length(unique(newtaomaxs))>1
        fprintf('maxtao detected\n');
        newarcthetas=repmat(0.02,1,cbnum);
    end
     
    cbr=cbr';
    cbr=cbr(:,1:end-1);
    
    
    [ ~,nlen,~,~] = adjdelta( cbr,newadcr,newdeltas,newtaomaxs,deltamin); 
    if len>=nlen
        len=nlen;
    end
    
    if isempty(len)==1
        len=deltalen;
    else
        deltalen=len;
    end
    
    arclength=arclength+len;
    
    datarec(Nxyz,arclength,newadcr,cbr,cbnum,newijkindex,btotal,xyz);
    
    cr=cbr;
    numcr=cbnum;
    deltas=newdeltas;
    adcr=newadcr;
    lost=newlost;
    ijkindex=newijkindex;
    arcthetas=newarcthetas;
    polanums=newpolanums;
    
    name=num2str(arclength);
    nameb=[path,name,hou];
    save(nameb,'cr','adcr','deltas','lost','ijkindex','arclength','arcthetas','polanums');
    save('arclen','arclength');
    savefig('acurrentbands.fig');
    
    fprintf('\narclength = %g\n ', arclength);
    [ deltas,len,con,numconvergence] = adjdelta( cr,adcr,deltas,newtaomaxs,deltamin); %adjust delta and detect convergence
    
    if numconvergence==1 %if convergence detected on all leaves, stop computation
        break;
    end
    
end
close(hbar);

end

