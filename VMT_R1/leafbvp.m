function [ br, taomax,delta,narctheta,polanum] = leafbvp( Nxyz,h,ti,i,numcr,cr,tao,adelta,deltao,adcr,deltamin,ijkindex,arcthetas,nullpoints,polanums,xyz,btotal )
%   boundry value problem for computing br
%   augments from left to right: function, timespan, current leaf, number 
%   of cr, all cr points, boundry value condition, arclength between cr 
%   and br, last cr point coresponds to current leaf
%   
% tao=300;
% adcri=adcr(i,:);

% deltamin=0.125;
% maxarc=0.7;

eps=1e-15;
odflag=0;
if length(adelta)==1
    odflag=1;
    delta=adelta;
else
    delta=adelta(i);
end
[ r,rl,rr,ijkcube ] = selectr( i,numcr,cr,ijkindex );  %select 3 points on Cr

polanum=zeros(3,1);
polanuminside=zeros(3,5);
% if r(1)>0&&r(1)<20&&r(2)>-7&&r(2)<7&&r(3)>-6.5&&r(3)<6.5
%     delta=delta/4;
% end

% deltao=0.5;
e=0.01;
if delta>deltao
    deltam=delta-deltao+deltao*(1+e);
    deltap=delta-deltao+deltao*(1-e);
else
    deltam=delta*(1+e);
    deltap=delta*(1-e);
end

stepnum=80;
% if numcr<=80
%     stepnum=20;
% elseif numcr>80&&numcr<=300
%     stepnum=40;
% elseif numcr>300&&numcr<=600
%     stepnum=40;
% elseif numcr>600
%     stepnum=50;
% end

stepmax=stepnum;
nde=norm(cr(:,i)'-adcr(i,:));
origin=[0,0,0];
dto=norm(cr(:,i)'-origin);


if nde<deltap||dto<4
    arctheta=arcthetas(i);
elseif delta<=deltamin
    arctheta=arcthetas(i);
%     arctheta=0.7;
else
    arctheta=0;
end

narctheta=arctheta;

nullflg=0;
cpnum=length(nullpoints(:,1));
for n=1:cpnum
    disrcp=norm(r-nullpoints(n,:));
    if disrcp<delta*1.5
        nullflg=1;
		break;
    end
end
if nullflg==1&&numcr>30
    arctheta=-0.9;
end


% arctheta4=arctheta;
oarctheta=arctheta;
% if arctheta>0.1
%    arctheta4=0; 
% end

if norm(r-origin)<2.8
    br=r;
    taomax=1;
    return;
end

nad= adnorm(i,numcr,adcr);
if norm(r-adcr(i,:))<deltamin*0.01
    br=r;
    taomax=1;
    return;
end
initao=11;
dirbr=adcr(i,:)-r; %growth direction vector for finding the next br
[ fn,ql,qr,qlijk,qrijk,qlstu,qrstu ] = geninit( Nxyz,r,rl,rr,initao,ijkcube,xyz ); % generate the init_points of 3 points
rstu=qrstu(1,:);
half=ceil(numcr);
idl=i-1;
idr=i+1;
extao=3;
oextao=extao;
mtao=11;
% taoc=((tao-1)/2)+1;
d=zeros(1,5);
dtaomax=d;
exbr=[r;r;r;r;r];
brtaomax=exbr;
taomax=0;
rv = asigncube( xyz, ijkcube );
[ ~,~,ijkcube ] = insidecell( Nxyz,r,rv,ijkcube,xyz );
[ti,h,con] = timestp( Nxyz,r,rstu,ti,h,ijkcube,delta,xyz,btotal ); % refine time steps

if con==1
    delta=0;
    br=r;
    taomax=1;
    return;
end

[index,tmpbr]=intersectfr( Nxyz,rstu,r,r,fn,ti,h,ijkcube,delta,arctheta,dirbr,xyz,btotal);
N=zeros(1,index);
deltas=repmat(delta,1,index);
for ij=1:index
    N(ij)=abs( norm(tmpbr(ij,:)-r)-deltas(ij));
end
[~,nn]=min(N);
d(1)=norm(r-tmpbr(nn,:));
dtaomax(1)=d(1);
exbr(1,:)=tmpbr(nn,:);
brtaomax(1,:)=exbr(1,:);

[flag,accbr] = accuracy( index,tmpbr,r,delta,dirbr,deltao,arctheta);
if flag==1
    br=accbr;
    return;
end

at=(tao-1)*(half-1)+1+(initao-1);
conl=zeros(1,at);
conr=zeros(1,at);
time=2;
while time<=half %generate all init_points around the cr
    stpl=0;
    stpr=0;
    il= anulus( idl,numcr);
    ir= anulus( idr,numcr);
    if odflag==0
        if adelta(il)==0
            stpl=1;
        end
        if adelta(ir)==0
            stpr=1;
        end
    end
    if stpl==1
        for ic=1:tao
            conl(length(ql(:,1))+ic-1)=1;
        end
    else
        if conl(length(ql(:,1))-tao+1)==1
            for ic=1:tao-1
                conl(length(ql(:,1))-ic+1)=0;
            end
        end
    end
    [ r2,rl2,rr2,ijkcube2 ] = selectr( il,numcr,cr,ijkindex );
    [ ~,qll,~,qllijk,~,qllstu,~ ] = geninit( Nxyz,r2,rl2,rr2,tao,ijkcube2,xyz );
    ql=cat(1,ql,qll(2:end,:));
    qlijk=cat(1,qlijk,qllijk(2:end,:));
    qlstu=cat(1,qlstu,qllstu(2:end,:));
    
    if stpr==1
        for ic=1:tao
            conr(length(qr(:,1))+ic-1)=1;
        end
    else
        if conr(length(qr(:,1))-tao+1)==1
            for ic=1:tao-1
                conr(length(qr(:,1))-ic+1)=0;
            end
        end
    end
    [ r22,rl22,rr22,ijkcube22 ] = selectr( ir,numcr,cr,ijkindex );
    [ ~,~,qrr,~,qrrijk,~,qrrstu ] = geninit( Nxyz,r22,rl22,rr22,tao,ijkcube22,xyz );
    qr=cat(1,qr,qrr(2:end,:));
    qrijk=cat(1,qrijk,qrrijk(2:end,:));
    qrstu=cat(1,qrstu,qrrstu(2:end,:));
    
    idl=idl-1;
    idr=idr+1;
    time=time+1;
end
taor=length(qr(:,1));
taol=length(ql(:,1));

% tao=length(ql(:,1));
% il= anulus( idl,numcr);
% ir= anulus( idr,numcr);
% if il~=ir
%     taoc=tao-taoc;
%     
%     ql=ql(1:taoc,:);
%     qlijk=qlijk(1:taoc,:);
%     qlstu=qlstu(1:taoc,:);
%     
%     qr=qr(1:taoc,:);
%     qrijk=qrijk(1:taoc,:);
%     qrstu=qrstu(1:taoc,:);
% end

[ ~,qln,qrn,qlnijk,qrnijk,qlnstu,qrnstu ] = geninit( Nxyz,r,ql(2,:),qr(2,:),300,ijkcube,xyz );

[indexl,tmpbrl]=intersectfr( Nxyz,qlnstu(2,:),qln(2,:),r,fn,ti,h,qlnijk(2,:),delta,arctheta,dirbr,xyz,btotal);
[indexr,tmpbrr]=intersectfr( Nxyz,qrnstu(2,:),qrn(2,:),r,fn,ti,h,qrnijk(2,:),delta,arctheta,dirbr,xyz,btotal);

chg=0;
tspan=ti(1,end);
tspan=-tspan;
bti=[0,tspan];

[bindexl,btmpbrl]=intersectfr( Nxyz,qlnstu(2,:),qln(2,:),r,fn,bti,h,qlnijk(2,:),delta,arctheta,dirbr,xyz,btotal);
[bindexr,btmpbrr]=intersectfr( Nxyz,qrnstu(2,:),qrn(2,:),r,fn,bti,h,qrnijk(2,:),delta,arctheta,dirbr,xyz,btotal);
  
if indexl>0
    tmpdirl=tmpbrl(1,:)-r;
    parcl=(dirbr*tmpdirl')/(norm(dirbr)*norm(tmpdirl));
    if parcl>0
        chg=1;
    end
    if indexr==0&&parcl>-0.5
        chg=1;
    end
    if bindexr>0
        tmpdirl=btmpbrr(1,:)-r;
        parcl=(dirbr*tmpdirl')/(norm(dirbr)*norm(tmpdirl));
        if parcl<0
            chg=1;
        end
    end
    
end
if indexr>0&&chg==0
    tmpdirr=tmpbrr(1,:)-r;
    parcr=(dirbr*tmpdirr')/(norm(dirbr)*norm(tmpdirr));
    if parcr>0
        chg=1;
    end
    if indexl==0&&parcr>-0.5
        chg=1;
    end
     if bindexl>0
        tmpdirl=btmpbrl(1,:)-r;
        parcl=(dirbr*tmpdirl')/(norm(dirbr)*norm(tmpdirl));
        if parcl<0
            chg=1;
        end
%         if bindexr==0&&parcl>-0.5
%             chg=1;
%         end
    end
end

% if chg==0
%     if r(1)>0&&r(1)<10&&r(2)>-6&&r(2)<6&&r(3)>-5.5&&r(3)<5.5
%         chg=1;
%     end
% end
% if chg==1
%     if r(1)>10
%         chg=2;
%     end
% end


[flagl,accbrl] = accuracy( indexl,tmpbrl,r,delta,dirbr,deltao,arctheta);
[flagr,accbrr] = accuracy( indexr,tmpbrr,r,delta,dirbr,deltao,arctheta);
if flagl==1
    br=accbrl;
    return;
end
if flagr==1
    br=accbrr;
    return;
end

closeil=0;
closeir=0;
optao=40;
while indexl==0&&indexr==0 % open area detected near the r

    [ ~,qll,qrr,qllijk,qrrijk,qllstu,qrrstu ] = geninit( Nxyz,r,ql(tao,:),qr(tao,:),optao,ijkcube,xyz );
    [indexl,tmpbrl]=intersectfr( Nxyz,qllstu(2,:),qll(2,:),r,fn,ti,h,qllijk(2,:),delta,arctheta,dirbr,xyz,btotal);
    [indexr,tmpbrr]=intersectfr( Nxyz,qrrstu(2,:),qrr(2,:),r,fn,ti,h,qrrijk(2,:),delta,arctheta,dirbr,xyz,btotal);
    ql=qll;
    qr=qrr;
    qlijk=qllijk;
    qrijk=qrrijk;
    qlstu=qllstu;
    qrstu=qrrstu;
    if indexl>0||indexr>0
        taol=length(ql(:,1));
        conl=zeros(1,taol);
        taor=length(qr(:,1));
        conr=zeros(1,taor);
        [flagl,accbrl] = accuracy( indexl,tmpbrl,r,delta,dirbr,deltao,arctheta);
        [flagr,accbrr] = accuracy( indexr,tmpbrr,r,delta,dirbr,deltao,arctheta);
        if flagl==1
            br=accbrl;
            return;
        end
        if flagr==1
            br=accbrr;
            return;
        end
        if indexl>0
            closeil=1;
        end
        if indexr>0
            closeir=1;
        end
        
        break;
    end
end
while indexl==0%&&index>1
    [ ~,ql,~,qlijk,~,qlstu,~ ] = geninit( Nxyz,r,ql(2,:),qr(2,:),optao,ijkcube,xyz );
    if norm(ql(end,:)-ql(1,:))<eps
       break; 
    end
    [indexl,tmpbrl]=intersectfr( Nxyz,qlstu(2,:),ql(2,:),r,fn,ti,h,qlijk(2,:),deltao,arctheta,dirbr,xyz,btotal);
    [flagl,accbrl] = accuracy( indexl,tmpbrl,r,delta,dirbr,deltao,arctheta);
    if flagl==1
        br=accbrl;
        return;
    end
    taol=length(ql(:,1));
    conl=zeros(1,taol);
    if indexl>0
        closeil=1;
    end
    if indexl==1
        if norm(tmpbrl-r)<eps*10
            indexl=0;
            break;
        end
    end
end
while indexr==0%&&index>1
    [ ~,~,qr,~,qrijk,~,qrstu ] = geninit( Nxyz,r,ql(2,:),qr(2,:),optao,ijkcube,xyz );
    if norm(qr(end,:)-qr(1,:))<eps
       break; 
    end
    [indexr,tmpbrr]=intersectfr( Nxyz,qrstu(2,:),qr(2,:),r,fn,ti,h,qrijk(2,:),deltao,arctheta,dirbr,xyz,btotal);
    [flagr,accbrr] = accuracy( indexr,tmpbrr,r,delta,dirbr,deltao,arctheta);
    if flagr==1
        br=accbrr;
        return;
    end
    taor=length(qr(:,1));
    conr=zeros(1,taor);
    if indexr>0
        closeir=1;
    end
    if indexr==1
        if norm(tmpbrr-r)<eps*10
            indexr=0;
            break;
        end
    end
end

% num2c=(tao-1)*8;
% lentaol=length(ql(:,1));
% if lentaol>num2c
%     lentaol=num2c;
% end
% lentaor=length(qr(:,1));
% if lentaor>num2c
%     lentaor=num2c;
% end
% 
% qls=ql(1:lentaol,:);
% qlijks=qlijk(1:lentaol,:);
% qlstus=qlstu(1:lentaol,:);
% conls=conl(1:lentaol);
% 
% qrs=qr(1:lentaor,:);
% qrijks=qrijk(1:lentaor,:);
% qrstus=qrstu(1:lentaor,:);
% conrs=conr(1:lentaor);

dflg=[0,tao];
extao2=6;
ml=polanums(1,i);%previous jishu in this r
mr=polanums(2,i);

% flagb=0;
if chg==1&&nullflg==0&&(ml>3||mr>3)
    [flagb,accbrb,d,exbr] = bkbvp(Nxyz,tao,index,r,fn,h,bti,d,exbr,ql,qr,extao2,qlijk,qrijk,qlstu,qrstu,conl,conr,deltamin,delta,deltao,dirbr,ijkcube,nad,xyz,btotal);
    if flagb==1
        br=accbrb;
        return;
    end
%     arctheta=0;
%     stepmax=50;
%     dflg(1)=1;
end

if dto<5&&dflg(1)==0
    dflg(1)=1;
end

minumflg=polanums(3,:);
pluflg=minumflg(i);
% minumflg;    %deleted flg of bands 1 delte 2 add

if indexl>0      %looking for br from r to rl
    brsl=-999;
    dli=d(2);
    bkl=0;
    precomflgl=0;
    
    if indexl==1&&ml>3&&ml<floor(numcr/2)&&closeil==0
        %bkl=1;
        ilp=i-ml;
        [lpplus,lpminus] = polarcheck( ilp,i,minumflg,numcr,tao );
        if pluflg<0
            lpplus=0;
            lpminus=0;
        end
        interjishul=initao+ml*(tao-1)+lpplus-lpminus;
        [brsl,bkl] = polarities( Nxyz,initao,dflg,interjishul,qlstu,ql,qlijk,r,fn,ti,h,delta,deltam,arctheta,dirbr,xyz,btotal );
        if any(conl(1:brsl)==1)
            interjishul2=0;
            for in=1:brsl
                if conl(in)==1&&(in+4)<=brsl
                    if conl(in+1)==1&&conl(in+2)==1&&conl(in+3)==1&&conl(in+4)==1
                        interjishul2=in-(tao-1);
                        break;
                    end
                end
            end
            if interjishul2>0
                 [brsl,bkl] = polarities( Nxyz,initao,dflg,interjishul2,qlstu,ql,qlijk,r,fn,ti,h,delta,deltam,arctheta,dirbr,xyz,btotal );
            end
        end
    end

    if precomflgl==0
        [ accbr,flagl,tmpbrl,jishul,bkl,anglel,brsl,arctheta,mbr] = exintersect(Nxyz,dflg,bkl,exbr(2,:),r,qlstu,ql,fn,ti,h,qlijk,taol,delta,deltao,dirbr,dli,arctheta,conl,nad,brsl,xyz,btotal);
       
        if chg==1&&jishul>(initao+(tao-1)*2)
            po=floor((jishul-initao)/(tao-1))+1;
            polanums(1,i)=po-1;
        else
            polanums(1,i)=0;
        end
        polanum=polanums(:,i);
        polanuminside(:,2)=polanum;
        if flagl==1
            if anglel>0.5
                narctheta=-0.9;
            end
            br=accbr;
            return;
        end
    end
    if (d(2)>deltam||d(2)<norm(tmpbrl-r))&&(dirbr*(tmpbrl-r)')/(norm(dirbr)*norm(tmpbrl-r))<arctheta
        d(2)=norm(tmpbrl-r);
        exbr(2,:)=tmpbrl;
        dd=abs(d(2)-delta);
        if dd<abs(dtaomax(2)-delta)
            dtaomax(2)=d(2);
            brtaomax(2,:)=exbr(2,:);
        end        
        
    end
    md=norm(mbr-r);
    if abs(md-delta)<abs(dtaomax(2)-delta)
        dtaomax(2)=md;
        brtaomax(2,:)=mbr;
    end
    
    if jishul>0
        if conl(jishul)==1
            stepmax=20;
            mtao=6;
        end
%         if ceil(length(ql(:,1))/10)<jishul
%             mtao=21;
%             oextao=6;
%         end
    end
    
    if flagl==2
        stepmax=10;
    elseif chg==1&&conr(tao)==1
        stepmax=10;
        mtao=21;
    end
    
    
    if jishul~=taol&&jishul>0 % open area detected
        extao=mtao;
        cons=zeros(1,extao);
        [ ~,qlm,~,qlmijk,~,qlmstu,~ ] = geninit( Nxyz,ql(jishul,:),ql(jishul+1,:),ql(1,:),extao,qlijk(jishul,:),xyz );
        cstep=0;
        while flagl==0&&cstep<stepmax
            [ accbr,flagl,tmpbrl,jishul2,bkl,anglel,brsl,arctheta,mbr] = exintersect( Nxyz,dflg,bkl,exbr(2,:),r,qlmstu,qlm,fn,ti,h,qlmijk,extao,delta,deltao,dirbr,d(2),arctheta,cons,nad,brsl,xyz,btotal);
            
            if anglel<-0.2&&nullflg==1&&jishul2>1
                if norm(tmpbrl-r)<eps
                    arctheta=narctheta;
                end
            end
            if flagl==1 %if br was found
                if anglel>0.5
                    narctheta=-0.9;
                end
                br=accbr;
                return;
            end
            if jishul2==0
                break;
            end
            if (d(2)>deltam||d(2)<norm(tmpbrl-r))&&(dirbr*(tmpbrl-r)')/(norm(dirbr)*norm(tmpbrl-r))<arctheta
                d(2)=norm(tmpbrl-r);
                exbr(2,:)=tmpbrl;
                dmaxflg=1;
                dd=abs(d(2)-delta);
                if dd<abs(dtaomax(2)-delta)
                    dtaomax(2)=d(2);
                    brtaomax(2,:)=exbr(2,:);
                end
            else
                dmaxflg=0;
            end
            
            md=norm(mbr-r);
            if abs(md-delta)<abs(dtaomax(2)-delta)
                dtaomax(2)=md;
                brtaomax(2,:)=mbr;
            end
            
            if (jishul2==0&&indexr>0)||jishul2==extao%||(bkl==0&&gapl<1e-8&&jishul2>1)
                break;
            end
            ql2=qlm;
            extao=oextao;
            cons=zeros(1,extao);
            [  ~,qlm,~,qlmijk,~,qlmstu,~ ] = geninit( Nxyz,ql2(jishul2,:),ql2(jishul2+1,:),ql2(1,:),extao,qlmijk(jishul2,:),xyz );
            if norm(qlm(1,:)-qlm(end,:))<eps&&dmaxflg==0&&cstep>40
                break;
            end
            cstep=cstep+1;
        end
    end
    
end

if arctheta~=oarctheta
    arctheta=oarctheta;
end

stepmax=stepnum;
mtao=11;
% oextao=3;

if indexr>0  %looking for br from r to rr
    brsr=-999;
    dri=d(3);
    bkr=0;
    precomflgr=0;
    
    if indexr==1&&mr>3&&mr<floor(numcr/2)&&closeir==0
        %         bkr=1;
        irp=i+mr;
        [rpplus,rpminus] = polarcheck( irp,i,minumflg,numcr,tao );
        if pluflg<0
            rpplus=0;
            rpminus=0;
        end
        interjishur=initao+mr*(tao-1)+rpplus-rpminus;
        [brsr,bkr] = polarities( Nxyz,initao,dflg,interjishur,qrstu,qr,qrijk,r,fn,ti,h,delta,deltam,arctheta,dirbr,xyz,btotal );
        if any(conr(1:brsr)==1)
            interjishur2=0;
            for in=1:brsr
                if conr(in)==1&&(in+4)<=brsr
                    if conr(in+1)==1&&conr(in+2)==1&&conr(in+3)==1&&conr(in+4)==1
                        interjishur2=in-(tao-1);
                        break;
                    end
                end
            end
            if interjishur2>0
                [brsr,bkr] = polarities( Nxyz,initao,dflg,interjishur2,qrstu,qr,qrijk,r,fn,ti,h,delta,deltam,arctheta,dirbr,xyz,btotal );
            end
        end
    end
    
    if precomflgr==0
        [ accbr,flagr,tmpbrr,jishur,bkr,angler,brsr,arctheta,mbr] = exintersect(Nxyz,dflg,bkr,exbr(3,:),r,qrstu,qr,fn,ti,h,qrijk,taor,delta,deltao,dirbr,dri,arctheta,conr,nad,brsr,xyz,btotal);
        
        if chg==1&&jishur>(initao+(tao-1)*2)
            po=floor((jishur-initao)/(tao-1))+1;
            polanums(2,i)=po-1;
        else
            polanums(2,i)=0;
        end
        polanum=polanums(:,i);
        polanuminside(:,3)=polanum;
        if flagr==1
            if angler>0.5
                narctheta=-0.9;
            end
            br=accbr;
            return;
        end
    end
    
    if (d(3)>deltam||d(3)<norm(tmpbrr-r))&&(dirbr*(tmpbrr-r)')/(norm(dirbr)*norm(tmpbrr-r))<arctheta
        d(3)=norm(tmpbrr-r);
        exbr(3,:)=tmpbrr;
        dd=abs(d(3)-delta);
        if dd<abs(dtaomax(3)-delta)
            dtaomax(3)=d(3);
            brtaomax(3,:)=exbr(3,:);
        end
    end
    
    md=norm(mbr-r);
    if abs(md-delta)<abs(dtaomax(3)-delta)
        dtaomax(3)=md;
        brtaomax(3,:)=mbr;
    end
    
%     if arctheta>oarctheta&&d(3)>deltamin*0.2&&d(3)<deltam
%         br= exbr(3,:);
%         delta=d(3);
%         narctheta=arctheta;
%         taomax=1;
%         return;
%     end
    
    if jishur>0
        if conr(jishur)==1
            stepmax=20;
            mtao=6;
        end
%         if jishur>tao*20
%             mtao=21;
%             oextao=6;
%         end
    end
    
    if flagr==2
        stepmax=10;
    elseif chg==1&&conl(tao)==1
        stepmax=10;
        mtao=21;
    end
    
    if jishur~=taor&&jishur>0 % open area detected
        extao=mtao;
        cons=zeros(1,extao);
        [ ~,~,qrm,~,qrmijk,~,qrmstu ] = geninit( Nxyz,qr(jishur,:),qr(1,:),qr(jishur+1,:),extao,qrijk(jishur,:),xyz );
        cstep=0;
        while flagr==0&&cstep<stepmax
            [ accbr,flagr,tmpbrr,jishur2,bkr,angler,brsr,arctheta,mbr] = exintersect(Nxyz,dflg,bkr,exbr(3,:),r,qrmstu,qrm,fn,ti,h,qrmijk,extao,delta,deltao,dirbr,d(3),arctheta,cons,nad,brsr,xyz,btotal);  
            
            if angler<-0.2&&nullflg==1&&jishur2>1
                if norm(tmpbrr-r)<eps
                    arctheta=narctheta;
                end
            end
            if flagr==1 %if br was found
                if angler>0.5
                    narctheta=-0.9;
                end
                br=accbr;
                return;
            end
            if jishur2==0
                break;
            end
            if (d(3)>deltam||d(3)<norm(tmpbrr-r))&&(dirbr*(tmpbrr-r)')/(norm(dirbr)*norm(tmpbrr-r))<arctheta
                d(3)=norm(tmpbrr-r);
                exbr(3,:)=tmpbrr;
                dmaxflg=1;
                dd=abs(d(3)-delta);
                if dd<abs(dtaomax(3)-delta)
                    dtaomax(3)=d(3);
                    brtaomax(3,:)=exbr(3,:);
                end
            else
                dmaxflg=0;
            end
            
            md=norm(mbr-r);
            if abs(md-delta)<abs(dtaomax(3)-delta)
                dtaomax(3)=md;
                brtaomax(3,:)=mbr;
            end
            
            if (jishur2==0&&indexl>0)||jishur2==extao%(bkr==0&&gapr<1e-8&&jishur2>1)||
                break;
            end
            qr2=qrm; 
            extao=oextao;
            cons=zeros(1,extao);
            [ ~,~,qrm,~,qrmijk,~,qrmstu ] = geninit( Nxyz,qr2(jishur2,:),qr2(1,:),qr2(jishur2+1,:),extao,qrmijk(jishur2,:),xyz );
            if norm(qrm(1,:)-qrm(end,:))<eps&&dmaxflg==0&&cstep>40
                break;
            end
            cstep=cstep+1;
        end
    end
    
    

end

% if chg==2&&numcr>30
%     [flagb,accbrb,~,~] = bkbvp(Nxyz,tao,index,r,fn,h,bti,d,exbr,ql,qr,extao2,qlijk,qrijk,qlstu,qrstu,conl,conr,deltamin,delta,deltao,dirbr,ijkcube,nad,xyz,btotal);
%     if flagb==1
%         br=accbrb;
%         return;
%     end
% end


dis=zeros(1,5);
mn=0;
for j=1:5
    tmpdir=brtaomax(j,:)-r;
    if norm(tmpdir)==0
        mn=mn+1;
        dis(mn)=norm(brtaomax(j,:)-r);
        continue;
    end
%     if (((dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir)))<arctheta)%&&(norm(exbr(j,:)-r)<deltap)
        mn=mn+1;
        dis(mn)=norm(brtaomax(j,:)-r);
%     else
%         mn=mn+1;
%         dis(mn)=0;
%     end
end
dis=dis(1:mn);
dis2=dis-delta;
for ii=1:5
   dis2(ii)=norm(dis2(ii)); 
end
[maxdis,~]=min(dis2);
br=brtaomax(dis2==maxdis,:);% convergence detected on this leaf, then br=br(tao_max);
br=br(1,:);

polanum=polanuminside(:,dis2==maxdis);
polanum=polanum(:,1);
% if norm(r-adcr(i,:))<deltap
%     br=r;
% end
% tmpdir=br-r;
% narctheta=((dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir)));
delta=norm(br-r);
if delta<deltap
    taomax=1;
end


end

