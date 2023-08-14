function [ newpolanums,new_cbr,new_cr,new_deltas,new_taomaxs,new_lost,new_ijkindex,newarcthetas] = adjcbr_new( Nxyz,br,parnum,cr,tao,ti,h, deltas,deltao,adcr,deltamin,taomaxs,lost,ijkindex,oijkindex,arcthetas,nullpoints,polanums,opolanums,con,xyz,btotal )
%   adjusting mesh points
%  

arcflg=1; %switch for the deltak admision, 1 on and 0 off.

cdeltas=deltas;
odeltas=deltas;
obrs=br;
% deltaf=max(unique(deltas));
deltafmin=deltao*0.5;           %minimum distance between adjacent meshpoints
deltafmax=deltao*1.5;           %maximum distance between adjacent meshpoints

nu=0;                    %number of mesh points after the processing
crt=cr';
numcr=length(cr(1,:));
polanums(3,:)=zeros(1,numcr);
bpolanums=polanums;
%proper deltak admission%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alphamin=0.2;
alphamax=0.3;
damin=0.1*deltao;
damax=1*deltao;
% smalls=zeros(1,numcr);
fprintf('admitting deltak...\n ');
addcr=adcr;
admflg=zeros(1,numcr);

for i=1:numcr
    if taomaxs(i)==1
        odeltas(i)=0;
    end
end
ods=odeltas(odeltas>=deltamin);
[minod,~]=min(ods);

admflgs=zeros(1,numcr);
parnum2=ceil(numcr/2);
if parnum==1
    parnum2=1;
end

if arcflg==1
    
    for iii=1:numcr
        if admflgs(iii)==1
            continue;
        end
        parn=parnum2+iii-1;
        if parn>numcr
            parn=numcr;
        end
        
        if parnum>1
            
            parfor in=iii:parn
                %     in=45;
                admflgs(in)=1;
                tmpbr=br(in,:);
                deltak=deltas(in);  % delta to be computed for the next circle
                cdelta=deltak;   %current delta to be accepted
                if lost(in)==1||(deltak<deltamin&&taomaxs(in)==1)
                    continue;
                end
                %     small=0;
                %     times=0;
                while 1
                    %         times=times+1;
                    if norm(crt(in,:)-tmpbr)<1e-4
                        break;
                    end
                    p=crt(in,:)+((crt(in,:)-tmpbr)/norm(crt(in,:)-tmpbr))*norm(crt(in,:)-adcr(in,:));
                    alphak=norm(p-adcr(in,:))/norm(crt(in,:)-adcr(in,:));
                    dak=deltak*alphak;
                    if alphak>=alphamax||dak>=damax
                        if deltak>deltamin
                            if (deltak/2)<=deltamin %minimum bound
                                deltak=deltamin;
                            else
                                deltak=deltak/2;
                            end
                            cdelta=deltak;
                            [tmpbr,~,~,~,polanum] = leafbvp( Nxyz,h,ti,in,numcr,cr,tao,deltak,deltao,addcr,deltamin,oijkindex,arcthetas,nullpoints,opolanums,xyz,btotal );
                            polanums(:,in)=polanum;
                            admflg(in)=1;
                            continue;
                        end
                    end
                    if alphak<alphamax&&dak<damax
                        
                        if alphak<alphamin&&dak<damin
                            deltak=deltak*2;
                        end
                        break;
                    end
                    if deltak<=deltamin
                        break;
                    end
                end
                
                cdeltas(in)=cdelta;
                deltas(in)=deltak;
                br(in,:)=tmpbr;
                
            end
        else
            for in=iii:parn
                %     in=45;
                admflgs(in)=1;
                tmpbr=br(in,:);
                deltak=deltas(in);  % delta to be computed for the next circle
                cdelta=deltak;   %current delta to be accepted
                if lost(in)==1||(deltak<deltamin&&taomaxs(in)==1)
                    continue;
                end
                %     small=0;
                %     times=0;
                while 1
                    %         times=times+1;
                    if norm(crt(in,:)-tmpbr)<1e-4
                        break;
                    end
                    p=crt(in,:)+((crt(in,:)-tmpbr)/norm(crt(in,:)-tmpbr))*norm(crt(in,:)-adcr(in,:));
                    alphak=norm(p-adcr(in,:))/norm(crt(in,:)-adcr(in,:));
                    dak=deltak*alphak;
                    if alphak>=alphamax||dak>=damax
                        if deltak>deltamin
                            if (deltak/2)<=deltamin %minimum bound
                                deltak=deltamin;
                            else
                                deltak=deltak/2;
                            end
                            cdelta=deltak;
                            [tmpbr,~,~,~,polanum] = leafbvp( Nxyz,h,ti,in,numcr,cr,tao,deltak,deltao,addcr,deltamin,oijkindex,arcthetas,nullpoints,opolanums,xyz,btotal );
                            polanums(:,in)=polanum;
                            admflg(in)=1;
                            continue;
                        end
                    end
                    if alphak<alphamax&&dak<damax
                        
                        if alphak<alphamin&&dak<damin
                            deltak=deltak*2;
                        end
                        break;
                    end
                    if deltak<=deltamin
                        break;
                    end
                end
                
                cdeltas(in)=cdelta;
                deltas(in)=deltak;
                br(in,:)=tmpbr;
                
            end
        end
        if parn==numcr
            break;
        end
        
    end
end

c2deltas=cdeltas;
for i=1:numcr
    if taomaxs(i)==1%&&c2deltas(i)<deltamin
        c2deltas(i)=0;
    end
end
delta1=c2deltas(c2deltas>=deltamin);
[mindelta,~]=min(delta1);
if isempty(delta1)==0

    for i=1:numcr
        if taomaxs(i)==1
            c2deltas(i)=mindelta;
        end
    end
    if length(unique(c2deltas))~=1&&any(admflg==1)
        
        for j=1:numcr
            if deltas(j)>=deltamin
                deltas(j)=mindelta;
            end
        end
        if minod>mindelta
            
            ensti=num2str(numcr);
            hbar2=waitbar(0,'recomputing..');
            fprintf('recomputing...  \n');
            adelta=deltas;
            bvpflg=zeros(1,numcr);
            for i=1:numcr
                sti=num2str(i);
                v=['recomputing...',sti,'/',ensti];
                waitbar(i/numcr,hbar2,v);
                
                if cdeltas(i)==deltas(i)||cdeltas(i)<deltamin||deltas(i)<deltamin||bvpflg(i)==1
                    continue;
                end
                parn=parnum+i-1;
                if parn>numcr
                    parn=numcr;
                end
                if parnum>1
                    parfor ii=i:parn %recompute cbr and all points can be accepted
                        % i=35;
                        if cdeltas(ii)==deltas(ii)||cdeltas(ii)<deltamin||deltas(ii)<deltamin||bvpflg(ii)==1
                            continue;
                        end
                        [accbr,taomax,~,~,polanum] = leafbvp( Nxyz,h,ti,ii,numcr,cr,tao,adelta,deltao,adcr,deltamin,oijkindex,arcthetas,nullpoints,opolanums,xyz,btotal ); %boundry value problem for computing br
                        br(ii,:)=accbr;
                        taomaxs(ii)=taomax;
                        polanums(:,ii)=polanum;
                        bvpflg(ii)=1;
                    end
                else
                    for ii=i:parn %recompute cbr and all points can be accepted
                        % i=35;
                        if cdeltas(ii)==deltas(ii)||cdeltas(ii)<deltamin||deltas(ii)<deltamin||bvpflg(ii)==1
                            continue;
                        end
                        [accbr,taomax,~,~,polanum] = leafbvp( Nxyz,h,ti,ii,numcr,cr,tao,adelta,deltao,adcr,deltamin,oijkindex,arcthetas,nullpoints,opolanums,xyz,btotal ); %boundry value problem for computing br
                        br(ii,:)=accbr;
                        taomaxs(ii)=taomax;
                        polanums(:,ii)=polanum;
                        bvpflg(ii)=1;
                    end
                end
                if parn==numcr
                    break;
                end
            end
            cdeltas=deltas;
            close(hbar2);
            if parnum>1
                for jj=1:numcr
                    rijk=oijkindex(jj,:);
                    rcube=asigncube( xyz,rijk );
                    [ ~,~,nijk ] = insidecell( Nxyz,br(jj,:),rcube,rijk,xyz );
                    ijkindex(jj,:)=nijk;
                end
            else
                for jj=1:numcr
                    rijk=oijkindex(jj,:);
                    rcube=asigncube( xyz,rijk );
                    [ ~,~,nijk ] = insidecell( Nxyz,br(jj,:),rcube,rijk,xyz );
                    ijkindex(jj,:)=nijk;
                end
            end
        end
    end
    
    if length(unique(c2deltas))==1&&any(taomaxs==1)
        for i=1:numcr
            if admflg(i)==1&&taomaxs(i)==1
                br(i)=obrs(i);
            end
        end
        polanums=bpolanums;
    end
    
    for i=1:numcr
        if taomaxs(i)==1
            deltas(i)=mindelta;
        end
    end
end


%%%%%%%mesh adaptation

newarcthetas=zeros(1,numcr);
newlost=newarcthetas;
cbr=zeros(numcr,3);
new_cr=cbr;
tmpdeltas=newarcthetas;
tcdeltas=tmpdeltas;
newtaomaxs=newarcthetas;
newijkindex=cbr;
indexadd=zeros(1,1);
dirindex=zeros(1,3);
newpolanums=zeros(3,numcr);

skipm = checkbr(cr,br,numcr,deltafmin);
skipm=[skipm,skipm(1)];
ijkindex=[ijkindex;ijkindex(1,:)];
arcthetas=[arcthetas,arcthetas(1)];
br=[br;br(1,:)];
numcr=numcr+1;
crt=[crt;crt(1,:)];
con=[con,con(1)];
% cr=[cr,cr(:,1)];
% adcr=[adcr;adcr(1,:)];
deltas=[deltas,deltas(1)];
cdeltas=[cdeltas,cdeltas(1)];
taomaxs=[taomaxs,taomaxs(1)];
lost=[lost,lost(1)];
ni=0;%number of mesh points to be added

polanums=[polanums,polanums(:,1)];


for i=2:numcr
%     i=2;
    if skipm(i-1)==1
        newarcthetas=newarcthetas(1:end-1);
        newlost=newlost(1:end-1);
        cbr=cbr(1:end-1,:);
        new_cr=new_cr(1:end-1,:);
        tmpdeltas=tmpdeltas(1:end-1);
        tcdeltas=tcdeltas(1:end-1);
        newtaomaxs=newtaomaxs(1:end-1);
        newijkindex=newijkindex(1:end-1,:);
        newpolanums=newpolanums(:,1:end-1);
        polanums(3,i)=polanums(3,i)+1;
        
        continue;
    end
    cflg=0;
    if con(i-1)==1&&con(i)==1
        cflg=1;
    end
    
    adjdis=abs(norm(br(i,:)-br(i-1,:)));
    if adjdis<deltafmin&&skipm(i-1)==0&&numcr>21&&cflg==0
        br(i,:)=br(i-1,:);
        crt(i,:)=crt(i-1,:);
        deltas(i)=deltas(i-1);
        cdeltas(i)=cdeltas(i-1);
        arcthetas(i)=arcthetas(i-1);
        taomaxs(i)=taomaxs(i-1);
        ijkindex(i,:)=ijkindex(i-1,:);
        skipm(i)=skipm(i-1);
        lost(i-1)=1;
        lost(i)=1;
        polanums(3,i-1)=polanums(3,i-1)+1;
        polanums(3,i)= polanums(3,i-1);
        
        newpolanums=newpolanums(:,1:end-1);
        newarcthetas=newarcthetas(1:end-1);
        newlost=newlost(1:end-1);
        cbr=cbr(1:end-1,:);
        new_cr=new_cr(1:end-1,:);
        tmpdeltas=tmpdeltas(1:end-1);
        tcdeltas=tcdeltas(1:end-1);
        newtaomaxs=newtaomaxs(1:end-1);
        newijkindex=newijkindex(1:end-1,:);

    elseif adjdis>deltafmax&&(skipm(i-1)==0||skipm(i-1)==2)&&cflg==0
        ni=ni+1;
        nu=nu+1;
        newarcthetas(nu)=arcthetas(i-1);
        newlost(nu)=lost(i-1);
        cbr(nu,:)=br(i-1,:);
        new_cr(nu,:)=crt(i-1,:);
        tmpdeltas(nu)=deltas(i-1);
        tcdeltas(nu)=cdeltas(i-1);
        newtaomaxs(nu)=taomaxs(i-1);
        newijkindex(nu,:)=ijkindex(i-1,:);
        newpolanums(:,nu)=polanums(:,i-1);%%%%%%
        
        rl=crt(i-1,:);
        rr=crt(i,:);
        centr=(rl+rr)/2;
        adrl=br(i-1,:);
        adrr=br(i,:);
        centadbr=(adrl+adrr)/2;
        dirbr=centr-centadbr;
        centpola=polanums(:,i-1);
        for ii=1:2
            if polanums(ii,i-1)<polanums(ii,i)&&polanums(ii,i-1)~=0
                centpola(ii)=polanums(ii,i-1);
            elseif polanums(ii,i-1)>polanums(ii,i)&&polanums(ii,i)~=0
                centpola(ii)=polanums(ii,i);
            else
                centpola(ii)=0;
            end
        end
        centpola(3)=-1;
        nu=nu+1;
        
        newpolanums=cat(2,newpolanums,[0;0;0]);
        newarcthetas=cat(2,newarcthetas,0);
        newlost=cat(2,newlost,0);
        cbr=cat(1,cbr,[0,0,0]);
        new_cr=cat(1,new_cr,[0,0,0]);
        tmpdeltas=cat(2,tmpdeltas,0);
        tcdeltas=cat(2,tcdeltas,0);
        newtaomaxs=cat(2,newtaomaxs,0);
        newijkindex=cat(1,newijkindex,[0,0,0]);
        if ni>1
            indexadd=cat(2,indexadd,0);
            dirindex=cat(1,dirindex,[0,0,0]);
        end
        
        indexadd(ni)=nu;%index of adding mesh position 
        dirindex(ni,:)=dirbr;
        cbr(nu,:)=centadbr;%zeros(1,3);
        newijkindex(nu,:)=ijkindex(i-1,:);
        new_cr(nu,:)=centr;
        newpolanums(:,nu)=centpola;
        if isempty(delta1)==0
            tmpdeltas(nu)=mindelta;
            tcdeltas(nu)=mindelta;
        else
            if deltas(i-1)~=deltas(i)
                if deltas(i-1)>deltas(i)
                    tmpdeltas(nu)=deltas(i-1);
                    tcdeltas(nu)=cdeltas(i-1);
                else
                    tmpdeltas(nu)=deltas(i);
                    tcdeltas(nu)=cdeltas(i);
                end
            else
                tmpdeltas(nu)=deltas(i-1);
                tcdeltas(nu)=cdeltas(i-1);
            end
        end
        newtaomaxs(nu)=0;
        newarcthetas(nu)=arcthetas(i-1);
        if lost(i-1)==1
            newlost(nu)=1;
        else
            newlost(nu)=0;
        end
    else
        nu=nu+1;
        newlost(nu)=0;
        cbr(nu,:)=br(i-1,:);
        new_cr(nu,:)=crt(i-1,:);
        tmpdeltas(nu)=deltas(i-1);
        newarcthetas(nu)=arcthetas(i-1);
        tcdeltas(nu)=cdeltas(i-1);
        newtaomaxs(nu)=taomaxs(i-1);
        newijkindex(nu,:)=ijkindex(i-1,:);
        newpolanums(:,nu)=polanums(:,i-1);
    end
end

if ni>0  %BVP computation for adding mesh points
    fprintf('adding mesh points...  ');
    addbrs=zeros(ni,3);
    taos=zeros(1,ni);
    rijks=zeros(ni,3);
    polas=zeros(3,ni);
    if parnum>1
        parfor j=1:ni
            id=indexadd(j);
            direc=dirindex(j,:);
            [centmp,centtaomax,rijk,polanum] = addmesh(Nxyz,tao,ti,h,new_cr,id,tcdeltas,deltao,direc,deltamin,newijkindex,newarcthetas,nullpoints,newpolanums,xyz,btotal);
            addbrs(j,:)=centmp;
            taos(j)=centtaomax;
            polas(:,j)=polanum;
            rcube=asigncube( xyz,rijk );
            [ ~,~,nijk ] = insidecell( Nxyz,centmp,rcube,rijk,xyz );
            rijks(j,:)=nijk;
        end
    else
        for j=1:ni
            id=indexadd(j);
            direc=dirindex(j,:);
            [centmp,centtaomax,rijk,polanum] = addmesh(Nxyz,tao,ti,h,new_cr,id,tcdeltas,deltao,direc,deltamin,newijkindex,newarcthetas,nullpoints,newpolanums,xyz,btotal);
            addbrs(j,:)=centmp;
            taos(j)=centtaomax;
            polas(:,j)=polanum;
            rcube=asigncube( xyz,rijk );
            [ ~,~,nijk ] = insidecell( Nxyz,centmp,rcube,rijk,xyz );
            rijks(j,:)=nijk;
        end
    end
    
    for j=1:ni
        id2=indexadd(j);
        if norm(new_cr(id2,:)-addbrs(j,:))<2e-14&&norm(new_cr(id2,:)-cbr(id2,:))>(tmpdeltas(id2)/2*(1+0.05))
            if norm(new_cr(id2,:)-cbr(id2,:))<(tmpdeltas(id2)*(1+0.05))
                continue;
            end
        end
        newpolanums(:,id2)=polas(:,j);
        cbr(id2,:)=addbrs(j,:);
        newijkindex(id2,:)=rijks(j,:);
        newtaomaxs(id2)=taos(j);
    end
end

new_cbr=[cbr;cbr(1,:)];
% cbnum=length(new_cbr);
new_deltas=tmpdeltas;
new_taomaxs=newtaomaxs;
new_lost=newlost;

nijkindex=newijkindex;
numcc=length(cbr(:,1));
for ii=1:numcc
stu = xyztostu( Nxyz,new_cbr(ii,1),new_cbr(ii,2),new_cbr(ii,3),newijkindex(ii,:),xyz );
[ ~,~,nijkindex(ii,:) ] = whichcell( stu,newijkindex(ii,:) );
end

new_ijkindex=nijkindex;

end

