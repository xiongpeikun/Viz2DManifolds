function [flag,br,d,exbr] = bkbvp( Nxyz,tao,index,r,fn,h,bti,d,exbr,ql,qr,extao,qlijk,qrijk,qlstu,qrstu,conl,conr,deltamin,delta,deltao,dirbr,ijkcube,nad,xyz,btotal)
%BKBVP Summary of this function goes here
%   Function of time-reversed orbit method

flag=0;
br=[];
cons=zeros(1,extao);

if delta<deltamin
    arctheta=-0.94;
elseif delta==deltamin
    arctheta=-0.94;
else
    arctheta=-0.94;
end

for i=1:8
    if conl(i)==1||conr(i)==1
        return;
    end
end
% if r(1)>0&&r(1)<10&&r(2)>-4&&r(2)<4&&r(3)>-6.5&&r(3)<6.5&&delta<=deltamin
%     arctheta=0.3;
% end

num2c=10+(tao-1)*3;
lentaol=length(ql(:,1));
if lentaol>num2c
    lentaol=num2c;
end
lentaor=length(qr(:,1));
if lentaor>num2c
    lentaor=num2c;
end
% lentao=ceil(lentao/5);
ql=ql(1:lentaol,:);
qlijk=qlijk(1:lentaol,:);
qlstu=qlstu(1:lentaol,:);
conl=conl(1:lentaol);

qr=qr(1:lentaor,:);
qrijk=qrijk(1:lentaor,:);
qrstu=qrstu(1:lentaor,:);
conr=conr(1:lentaor);

[ ~,qln,qrn,qlnijk,qrnijk,qlnstu,qrnstu ] = geninit( Nxyz,r,ql(2,:),qr(2,:),100,ijkcube,xyz );
[indexl,tmpbrl]=intersectfr( Nxyz,qlnstu(2,:),qln(2,:),r,fn,bti,h,qlnijk(2,:),delta,arctheta,dirbr,xyz,btotal);
[indexr,tmpbrr]=intersectfr( Nxyz,qrnstu(2,:),qrn(2,:),r,fn,bti,h,qrnijk(2,:),delta,arctheta,dirbr,xyz,btotal);
[flagl,accbrl] = accuracy( indexl,tmpbrl,r,delta,dirbr,deltao,arctheta);
[flagr,accbrr] = accuracy( indexr,tmpbrr,r,delta,dirbr,deltao,arctheta);
if flagl==1
    br=accbrl;
    flag=1;
    return;
end
if flagr==1
    br=accbrr;
    flag=1;
    return;
end

taol=length(ql(:,1));
taor=length(qr(:,1));

optao=40;
while indexl==0&&indexr==0 % open area detected near the r

    [ ~,qll,qrr,qllijk,qrrijk,qllstu,qrrstu ] = geninit( Nxyz,r,ql(2,:),qr(2,:),optao,ijkcube,xyz );
    [indexl,tmpbrl]=intersectfr( Nxyz,qllstu(2,:),qll(2,:),r,fn,bti,h,qllijk(2,:),delta,arctheta,dirbr,xyz,btotal);
    [indexr,tmpbrr]=intersectfr( Nxyz,qrrstu(2,:),qrr(2,:),r,fn,bti,h,qrrijk(2,:),delta,arctheta,dirbr,xyz,btotal);
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
            flag=1;
            return;
        end
        if flagr==1
            br=accbrr;
            flag=1;
            return;
        end
        break;
    end
end
while indexl==0&&index>1
    [ ~,ql,~,qlijk,~,qlstu,~ ] = geninit( Nxyz,r,ql(2,:),qr(2,:),optao,ijkcube,xyz );
    [indexl,tmpbrl]=intersectfr( Nxyz,qlstu(2,:),ql(2,:),r,fn,bti,h,qlijk(2,:),delta,arctheta,dirbr,xyz,btotal);
    [flagl,accbrl] = accuracy( indexl,tmpbrl,r,delta,dirbr,deltao,arctheta);
    if flagl==1
        flag=1;
        br=accbrl;
        return;
    end
    taol=length(ql(:,1));
    conl=zeros(1,taol);
end
while indexr==0&&index>1
    [ ~,~,qr,~,qrijk,~,qrstu ] = geninit( Nxyz,r,ql(2,:),qr(2,:),optao,ijkcube,xyz );
    [indexr,tmpbrr]=intersectfr( Nxyz,qrstu(2,:),qr(2,:),r,fn,bti,h,qrijk(2,:),delta,arctheta,dirbr,xyz,btotal);
    [flagr,accbrr] = accuracy( indexr,tmpbrr,r,delta,dirbr,deltao,arctheta);
    if flagr==1
        flag=1;
        br=accbrr;
        return;
    end
    taor=length(qr(:,1));
    conr=zeros(1,taor);
end
dflg=[0,tao];
  
if indexl>0      %looking for br from r to rl
    
    dli=d(4);
    bkl=0;
    brsl=-999;
    [ accbr,flagl,tmpbrl,bjishul,bkl,arcthetal,brsl,arctheta] = exintersect(Nxyz,dflg,bkl,exbr(4,:),r,qlstu,ql,fn,bti,h,qlijk,taol,delta,deltao,dirbr,dli,arctheta,conl,nad,brsl,xyz,btotal);
    
    if flagl==1
        flag=1;
        br=accbr;
        return;
    end
    
    exbr(4,:)=tmpbrl;
    d(4)=norm(tmpbrl-r);
    if isnan(arcthetal)
        arcthetal=arctheta;
    else
        arcthetal=arctheta;
    end
    
    if bjishul~=taol&&bjishul>0&&conl(bjishul)==0 % open area detected
        [ ~,qlm,~,qlmijk,~,qlmstu,~ ] = geninit( Nxyz,ql(bjishul,:),ql(bjishul+1,:),ql(1,:),extao,qlijk(bjishul,:),xyz );
        cstep=0;
        while flagl==0&&cstep<100
            [ accbr,flagl,tmpbrl,bjishul2,bkl,~,brsl] = exintersect( Nxyz,dflg,bkl,exbr(4,:),r,qlmstu,qlm,fn,bti,h,qlmijk,extao,delta,deltao,dirbr,d(4),arcthetal,cons,nad,brsl,xyz,btotal);
            
            if flagl==1 %if br was found
                flag=1;
                br=accbr;
                return;
            end
            if bjishul2==0
                break;
            end
            if d(4)<norm(tmpbrl-r)
                d(4)=norm(tmpbrl-r);
                exbr(4,:)=tmpbrl;
            end
            if (bjishul2==0&&indexr>0)||bjishul2==extao
                break;
            end
            ql2=qlm;
            [  ~,qlm,~,qlmijk,~,qlmstu,~ ] = geninit( Nxyz,ql2(bjishul2,:),ql2(bjishul2+1,:),ql2(1,:),extao,qlmijk(bjishul2,:),xyz );
            cstep=cstep+1;
        end
    end
    
end

if indexr>0  %looking for br from r to rr
    
    dri=d(5);
    bkr=0;
    brsr=-999;
    [ accbr,flagr,tmpbrr,bjishur,bkr,arcthetar,brsr,arctheta] = exintersect(Nxyz,dflg,bkr,exbr(5,:),r,qrstu,qr,fn,bti,h,qrijk,taor,delta,deltao,dirbr,dri,arctheta,conr,nad,brsr,xyz,btotal);
    
    if flagr==1
        flag=1;
        br=accbr;
        return;
    end
    
    exbr(5,:)=tmpbrr;
    d(5)=norm(tmpbrr-r);
    if isnan(arcthetar)
        arcthetar=arctheta;
    else
        arcthetar=arctheta;
    end
    
    
    if bjishur~=taor&&bjishur>0&&conr(bjishur)==0 % open area detected
        [ ~,~,qrm,~,qrmijk,~,qrmstu ] = geninit( Nxyz,qr(bjishur,:),qr(1,:),qr(bjishur+1,:),extao,qrijk(bjishur,:),xyz );
        cstep=0;
        while flagr==0&&cstep<100
            [ accbr,flagr,tmpbrr,bjishur2,bkr,~,brsr] = exintersect(Nxyz,dflg,bkr,exbr(5,:),r,qrmstu,qrm,fn,bti,h,qrmijk,extao,delta,deltao,dirbr,d(5),arcthetar,cons,nad,brsr,xyz,btotal);
            
            if flagr==1 %if br was found
                flag=1;
                br=accbr;
                return;
            end
            if bjishur2==0
                break;
            end
            
            if d(5)<norm(tmpbrr-r)
                d(5)=norm(tmpbrr-r);
                exbr(5,:)=tmpbrr;
            end
            if (bjishur2==0&&indexl>0)||bjishur2==extao
                break;
            end
            qr2=qrm;
            [ ~,~,qrm,~,qrmijk,~,qrmstu ] = geninit( Nxyz,qr2(bjishur2,:),qr2(1,:),qr2(bjishur2+1,:),extao,qrmijk(bjishur2,:),xyz );
            cstep=cstep+1;
        end
    end
    
end





end

