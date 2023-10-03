function [ br,flag,exbr,jishu,bk,arctheta3,brs,arctheta,mbr] = exintersect( Nxyz,dflg2,bk,exbr,r,qmstu,qm,fn,ti,h,qmijk,extao,delta,deltao,dirbr,d1,arctheta,con,nad,tmpbr,xyz,btotal)
%EXINTERSECT Summary of this function goes here
%   Function of Bisection Method

dflg=dflg2(1);
tao=dflg2(2);
arctheta2=0.99;

e=0.01;
eps=2e-14;
if delta>deltao
    deltam=delta-deltao+deltao*(1+e);
    deltap=delta-deltao+deltao*(1-e);
else
    deltam=delta*(1+e);
    deltap=delta*(1-e);
end
md=0;

d2=0;
flag=0;
br=[];
jishu=0;
di=d1;
nn=1;
bnn=nn;

brs=exbr;
fircon=0;

arctheta3=arctheta;
okflg=0;
if length(nad(1,:))==3
    okflg=1;
end
conjishu=0;
dis=0;
disbrs=r;
skipflg=0;
if length(tmpbr(1,:))==3
    pts=length(tmpbr(:,1));
elseif length(tmpbr)==1&&tmpbr>0
    i=tmpbr;
    skipflg=1;
end
arcchg=0;
connum=(tao-1)*2;
mbr=r;

for m=1:extao
    bdis=dis;
    bdbrs=disbrs;
    larctheta3=arctheta3;
    
    if con(m)==1&&fircon==0
        conjishu=m;
        fircon=1;
    end
    if fircon==1
        arctheta=-0.7;
    end
    if con(m)==1&&m>1&&connum==0
        jishu=m-1;
        nn=bnn;
        bk=bbk;
        exbr=brs(nn,:);
        arctheta3=arctheta;
        return;
    end
    if con(m)==1&&connum>0
        connum=connum-1;
    end
    if connum>0&&connum<tao*2&&con(m)==0
        connum=(tao-1)*2;
    end
    if fircon==2
        fircon=0;
    end
    if fircon==1&&conjishu~=m
        fircon=2;
    end
    
    bbk=bk;
    if m==1&&length(tmpbr(1,:))==3
        
    elseif skipflg==1&&m>1&&m<=i
        continue;
    else
        [pts,tmpbr,oflg]=intersectfr( Nxyz,qmstu(m,:),qm(m,:),r,fn,ti,h,qmijk(m,:),delta,arctheta,dirbr,xyz,btotal);
    end
    
    if pts==0
        if oflg==1&&dflg==1&&bbk==1
            if fircon==2
                 jishu=conjishu;
                 return;
            end
            continue;
        else
            if fircon==2
                jishu=conjishu;
            else
                jishu=m-1;
            end
            nn=bnn;
            bk=bbk;
            arctheta3=larctheta3;
            if abs(md-delta)<abs(norm(brs(nn,:)-r)-delta)
                exbr=mbr;
            else
                exbr=brs(nn,:);
            end
            if bbk==1
               exbr=r; 
            end
            return;
        end
    end

    chgflg=0;
    dst=zeros(1,pts);
    for ij=1:pts
        tmpdir=tmpbr(ij,:)-r;
        if ((dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir))<arctheta2)
            dst(ij)=norm(tmpdir);
        elseif ((dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir))>=arctheta2)&&norm(tmpdir)>deltam
            dst(ij)=0;
            tmpbr(ij,:)=r;
        end
    end
    deltas=repmat(delta,1,pts);
    lastd=repmat(di,1,pts);
    [~,nn1]=min(abs(dst-deltas));
    [~,nn]=min(abs(dst-lastd));
    cg=0;
    if nn~=nn1
        nn=nn1;
        chgflg=1;
        cg=1;
    end
    if nn==nn1&&length(brs(:,1))~=pts&&cg==0
        if ((length(brs(:,1))-bnn)==(pts-nn)&&chgflg==0&&nn>bnn)||nn==bnn
            chgflg=0;
        else
            chgflg=1;
        end
    end
    if bnn~=nn&&length(brs(:,1))==pts
        chgflg=1;
    end

    d2=norm(tmpbr(nn,:)-r);
    dis=d2;
    disbrs=tmpbr(nn,:);
    tmpdir=tmpbr(nn,:)-r;
    arctheta3=(dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir));
    
    if bbk==1&&(arctheta3-larctheta3)>-0.15&&arctheta3<arctheta&&d2>delta*2
        brs=tmpbr;
        di=d2;
        bnn=nn;
        continue;
    end
    
    if ((isnan(larctheta3)==1||larctheta3>0.8)&&d2>deltam*2&&arctheta3<0.01)&&m>1&&arctheta3<0.7
        if isnan(larctheta3)==1||(arctheta3-larctheta3)<-0.7
            [ ~,~,qpm,~,qpmijk,~,qpmstu ] = geninit( Nxyz,qm(m-1,:),qm(1,:),qm(m,:),41,qmijk(m,:),xyz );
            [ptspb,tmpbrpb]=intersectfr( Nxyz,qpmstu(38,:),qpm(38,:),r,fn,ti,h,qpmijk(38,:),delta,arctheta,dirbr,xyz,btotal);
            [flag,accbr] = accuracy( ptspb,tmpbrpb,r,delta,dirbr,deltao,arctheta );
            if flag==1
                br=accbr;
                jishu=m-1;
                return;
            end
            dstp=zeros(1,ptspb);
            for ij=1:ptspb
                tmpdir=tmpbrpb(ij,:)-r;
                if ((dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir))<arctheta2)
                    dstp(ij)=norm(tmpdir);
                elseif ((dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir))>=arctheta2)
                    dstp(ij)=0;
                    tmpbrpb(ij,:)=r;
                end
            end
            lastd=repmat(d2,1,ptspb);
            [~,nnpb]=min(abs(dstp-lastd));
            dnum=0;
            while (ptspb~=pts||nnpb~=nn)&&dnum<10
                [ ~,~,qpm,~,qpmijk,~,qpmstu ] = geninit( Nxyz,qpm(9,:),qpm(1,:),qpm(10,:),11,qpmijk(1,:),xyz );
                [ptspb,tmpbrpb]=intersectfr( Nxyz,qpmstu(5,:),qpm(5,:),r,fn,ti,h,qpmijk(5,:),delta,arctheta,dirbr,xyz,btotal);
                dstp=zeros(1,ptspb);
                for ij=1:ptspb
                    tmpdir=tmpbrpb(ij,:)-r;
                    if ((dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir))<arctheta2)
                        dstp(ij)=norm(tmpdir);
                    elseif ((dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir))>=arctheta2)
                        dstp(ij)=0;
                        tmpbrpb(ij,:)=r;
                    end
                end
                lastd=repmat(d2,1,ptspb);
                [~,nnpb]=min(abs(dstp-lastd));
                dnum=dnum+1;
            end
            if ptspb==0%%%%%%
                jishu=0;
                return;
            end
            if ptspb==pts
                d2pb=norm(tmpbrpb(nnpb,:)-r);
                crementpb=d2-d2pb;
                
                tmpdirp=tmpbrpb(nn,:)-r;
                arctheta3pb=(dirbr*tmpdirp')/(norm(dirbr)*norm(tmpdirp));
                
                if arctheta3<arctheta&&arctheta3pb<0&&d2>deltam
                    jishu=m-1;
                    arctheta3=larctheta3;
                    nn=bnn;
                    bk=bbk;
                    exbr=brs(nn,:);
                    return;
                end
                
                if arctheta3pb<arctheta3&&crementpb>0
                    arctheta3=arctheta-0.001;
                end
            end
        end
    end
    
    if pts~=0
        [flag,accbr] = accuracy( pts,tmpbr,r,delta,dirbr,deltao,arctheta );
        if flag==1
            jishu=m-1;
            br=accbr;
            return;
        end
    end
    
    if okflg==1&&pts==1&&arctheta3>0&&arcchg==0%%%%2021211
        adcr=dirbr+r;
        ndir1=tmpbr-adcr;
        ndir2=tmpbr-nad(5,:);
        ndir3=tmpbr-nad(6,:);
        jia=norm(dot(ndir1,nad(1,:)));
        jia2=norm(dot(ndir2,nad(2,:)));
        jia3=norm(dot(ndir3,nad(3,:)));
        
        if (jia<0.0001||jia2<0.0001||jia3<0.0001)&&dis<deltam&&m>1%&&jia2<0.5
            if bbk==1
                brs=r;
                bnn=1;
            else
                brs=tmpbr;
                bnn=nn;
            end
            if (arctheta3<arctheta)
                bk=0;
            else
                bk=1;
            end
            continue;
        end
        if (dis>deltam&&bdis>deltam&&dis-bdis>0)||(dis>deltam*2&&bdis>deltam)
            if bbk==1
                brs=r;
                di=0;
                bnn=1;
            else
                brs=tmpbr;
                di=d2;
                bnn=nn;
            end
            if (arctheta3<arctheta)
                bk=0;
            else
                bk=1;
            end
            continue;
        end
    end
    if (arctheta3<arctheta)||(m==1&&d2<eps)
        if bbk==1&&arctheta3>0
            if ((dis>deltam&&bdis>deltam)||(dis>deltam*2&&bdis>deltam))&&pts==1
                if (((bdbrs-r)*tmpdir')/(norm(bdbrs-r)*norm(tmpdir))>0)
                    brs=r;
                    di=0;
                    bnn=1;
                    continue;
                end
            end
            di=0;
        end
        bk=0;
        minter=abs(d2-delta);
        mintero=abs(md-delta);
        crement=d2-di;
        if minter<mintero
            if crement>0&&d2<deltap||crement<0&&d2>deltam
                mbr=tmpbr(nn,:);
                md=norm(mbr-r);
            end
        end
        exbr=tmpbr(nn,:);
        
        if chgflg==1&&length(brs(:,1))>1&&pts==1&&d2>deltam
            jishu=m-1;
            nn=bnn;
            bk=bbk;
            exbr=brs(nn,:);
            arctheta3=larctheta3;
            return;
        end
        if (chgflg==1&&m==extao&&cg==0)
            jishu=m-1;
            nn=bnn;
            bk=bbk;
            arctheta3=larctheta3;
            if abs(md-delta)<abs(norm(brs(nn,:)-r)-delta)
                exbr=mbr;
            else
                exbr=brs(nn,:);
            end
            return;
        end
		
        if chgflg==1&&cg==1&&m==extao
            [ ~,~,qpm,~,qpmijk,~,qpmstu ] = geninit( Nxyz,qm(m-1,:),qm(1,:),qm(m,:),10,qmijk(m,:),xyz );
            [ptspb,tmpbrpb]=intersectfr( Nxyz,qpmstu(7,:),qpm(7,:),r,fn,ti,h,qpmijk(7,:),delta,arctheta,dirbr,xyz,btotal);
            [flag,accbr] = accuracy( ptspb,tmpbrpb,r,delta,dirbr,deltao,arctheta );
            if flag==1
                br=accbr;
                jishu=m-1;
                return;
            end
            dstp=zeros(1,ptspb);
            for ij=1:ptspb
                tmpdir=tmpbrpb(ij,:)-r;
                if ((dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir))<arctheta2)
                    dstp(ij)=norm(tmpdir);
                elseif ((dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir))>=arctheta2)
                    dstp(ij)=0;
                    tmpbrpb(ij,:)=r;
                end
            end
            lastd=repmat(d2,1,ptspb);
            [~,nnpb]=min(abs(dstp-lastd));
            dnum=0;
            while (ptspb~=pts||nnpb~=nn)&&dnum<10
                [ ~,~,qpm,~,qpmijk,~,qpmstu ] = geninit( Nxyz,qpm(9,:),qpm(1,:),qpm(10,:),11,qpmijk(1,:),xyz );
                [ptspb,tmpbrpb]=intersectfr( Nxyz,qpmstu(5,:),qpm(5,:),r,fn,ti,h,qpmijk(5,:),delta,arctheta,dirbr,xyz,btotal);
                dstp=zeros(1,ptspb);
                for ij=1:ptspb
                    tmpdir=tmpbrpb(ij,:)-r;
                    if ((dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir))<arctheta2)
                        dstp(ij)=norm(tmpdir);
                    elseif ((dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir))>=arctheta2)
                        dstp(ij)=0;
                        tmpbrpb(ij,:)=r;
                    end
                end
                lastd=repmat(d2,1,ptspb);
                [~,nnpb]=min(abs(dstp-lastd));
                dnum=dnum+1;
            end
            if ptspb==0%%%%%%
                jishu=0;
                return;
            end
            d2pb=norm(tmpbrpb(nnpb,:)-r);
            crementpb=d2-d2pb;
            if any(dstp>0)
                if (d2<deltap&&crementpb<eps)||(crementpb>-eps&&d2>deltam)
                    jishu=m-1;
                    arctheta3=larctheta3;
                    nn=bnn;
                    bk=bbk;
                    exbr=brs(nn,:);
                    return;
                end
            end
            
        end
        
        if ((chgflg==1&&m<extao)||(m==1&&pts>1))
            
            [ ~,~,qpm,~,qpmijk,~,qpmstu ] = geninit( Nxyz,qm(m,:),qm(1,:),qm(m+1,:),11,qmijk(m,:),xyz );
            nannum=length(find(isnan(qpm(:,1))==1));%%%%%%jingduxiaxian
            if nannum~=0
                err=norm((d2-delta)/delta);
                if err<0.05
                    br=tmpbr(nn,:);
                    flag=1;
                    return;
                end
            end
            [ptsp,tmpbrp]=intersectfr( Nxyz,qpmstu(5,:),qpm(5,:),r,fn,ti,h,qpmijk(5,:),delta,arctheta,dirbr,xyz,btotal);
            dstp=zeros(1,ptsp);
            for ij=1:ptsp
                tmpdir=tmpbrp(ij,:)-r;
                if ((dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir))<arctheta2)
                    dstp(ij)=norm(tmpdir);
                elseif ((dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir))>=arctheta2)
                    dstp(ij)=0;
                    tmpbrp(ij,:)=r;
                end
            end
            lastd=repmat(d2,1,ptsp);
            [~,nnp1]=min(abs(dstp-lastd));
            dnum=0;
            while (ptsp~=pts||nnp1~=nn)&&dnum<20
                [ ~,~,qpm,~,qpmijk,~,qpmstu ] = geninit( Nxyz,qpm(1,:),qpm(1,:),qpm(2,:),11,qpmijk(1,:),xyz );
                [ptsp,tmpbrp]=intersectfr( Nxyz,qpmstu(5,:),qpm(5,:),r,fn,ti,h,qpmijk(5,:),delta,arctheta,dirbr,xyz,btotal);
                dstp=zeros(1,ptsp);
                for ij=1:ptsp
                    tmpdir=tmpbrp(ij,:)-r;
                    if ((dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir))<arctheta2)
                        dstp(ij)=norm(tmpdir);
                    elseif ((dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir))>=arctheta2)
                        dstp(ij)=0;
                        tmpbrp(ij,:)=r;
                    end
                end
                lastd=repmat(d2,1,ptsp);
                [~,nnp1]=min(abs(dstp-lastd));
                dnum=dnum+1;
            end
            if isempty(nnp1)==1
                jishu=0;
                br=exbr;
                return;
            end
            d2p=norm(tmpbrp(nnp1,:)-r);
            crementp=d2p-d2;
            
            if bbk==1&&m>1&&any(dstp>0)
                if d2>deltam&&crementp>-eps
                    jishu=m-1;
                    arctheta3=larctheta3;
                    bk=bbk;
                    return;
                end
            end
            
            if m<extao&&any(dstp>0)
                if crement>-eps&&crementp>-eps&&d2<deltap||crement<eps&&crementp<eps&&d2>deltam
                    brs=tmpbr;
                    di=d2;
                    bnn=nn;
                    continue;
                end
                if nn~=bnn
                    if crementp>-eps&&d2<deltap||crementp<eps&&d2>deltam
                        brs=tmpbr;
                        di=d2;
                        bnn=nn;
                        continue;
                    end
                end
            end
            
            if m>1&&crementp>-eps&&(d2>deltam||crement<eps&&d2<deltap)&&any(dstp>0)
                jishu=m-1;
                arctheta3=larctheta3;
                bk=bbk;
                return;
            end
            
            if m>1&&any(dstp>0)
                crementpp=d2p-norm(brs(bnn,:)-r);%%%%%%%%%%%%%%%%%%
                if (crementpp<eps&&d2<deltap&&d2p>deltam)||(crementpp>-eps&&d2>deltam&&d2p<deltap)
                    jishu=m;
                    brs=tmpbr;
                    return;
                end
            end
            
            if (d2<deltap&&d2p>deltam||d2>deltam&&d2p<deltap)&&any(dstp>0)
                jishu=m;
                brs=tmpbr;
                return;
            end
            
            [ptsp2,tmpbrp2]=intersectfr( Nxyz,qmstu(m+1,:),qm(m+1,:),r,fn,ti,h,qmijk(m+1,:),delta,arctheta,dirbr,xyz,btotal);
            dstp=zeros(1,ptsp2);
            for ij=1:ptsp2
                tmpdir=tmpbrp2(ij,:)-r;
                if ((dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir))<arctheta2)
                    dstp(ij)=norm(tmpdir);
                elseif ((dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir))>=arctheta2)
                    dstp(ij)=0;
                    tmpbrp2(ij,:)=r;
                end
            end
            lastd=repmat(d2,1,ptsp2);
            deltas=repmat(delta,1,ptsp2);
            [~,nnp]=min(abs(dstp-lastd));
            [~,nnp2]=min(abs(dstp-deltas));
            if isempty(nnp)==1
                if m==1
                    jishu=m;
                    brs=tmpbr;
                    return;
                end

                if ((crementp<eps&&d2>deltam&&d2p>deltam)||(crementp>-eps&&d2<deltap&&d2p<deltap))
                    brs=tmpbr;
                    di=d2;
                    bnn=nn;
                    continue;
                end
                if (crementp<eps&&d2>deltam&&d2p<deltap)||(crementp>-eps&&d2<deltap&&d2p>deltam)
                    jishu=m;
                    brs=tmpbr;
                    return;
                end
                [flag,accbr] = accuracy( ptsp,tmpbrp,r,delta,dirbr,deltao,arctheta );
                if flag==1
                    br=accbr;
                    jishu=m-1;
                    return;
                end
                
                jishu=m-1;
                nn=bnn;
                bk=bbk;
                exbr=brs(nn,:);
                return;
            end
            if isempty(nnp)==0&&any(dstp>0)
                [flag,accbr] = accuracy( ptsp2,tmpbrp2,r,delta,dirbr,deltao,arctheta );
                if flag==1
                    br=accbr;
                    jishu=m-1;
                    return;
                end
                cg=0;
                if nnp2~=nnp
                    tnn=nnp;
                    nnp=nnp2;
                    nnp2=tnn;
                    cg=1;
                end
                d2p2=norm(tmpbrp2(nnp,:)-r);
                crementp2=d2p2-d2p;
                
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
               if m<extao-1
                   [ ~,~,qpm,~,qpmijk,~,qpmstu ] = geninit( Nxyz,qm(m+1,:),qm(1,:),qm(m+2,:),11,qmijk(m,:),xyz );
                   nannum=length(find(isnan(qpm(:,1))==1));%%%%%%jingduxiaxian
                   if nannum~=0
                       err=norm((d2-delta)/delta);
                       if err<0.05
                           br=tmpbr(nn,:);
                           flag=1;
                           return;
                       end
                   end
                   [ptsp3,tmpbrp3]=intersectfr( Nxyz,qpmstu(5,:),qpm(5,:),r,fn,ti,h,qpmijk(5,:),delta,arctheta,dirbr,xyz,btotal);
                   dstp=zeros(1,ptsp3);
                   for ij=1:ptsp3
                       tmpdir=tmpbrp3(ij,:)-r;
                       if ((dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir))<arctheta2)
                           dstp(ij)=norm(tmpdir);
                       elseif ((dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir))>=arctheta2)
                           dstp(ij)=0;
                           tmpbrp3(ij,:)=r;
                       end
                   end
                   lastd=repmat(d2,1,ptsp3);
                   [~,nnp3]=min(abs(dstp-lastd));
                   dnum=0;
                   while (ptsp3~=ptsp2||nnp3~=nnp2)&&dnum<20
                       [ ~,~,qpm,~,qpmijk,~,qpmstu ] = geninit( Nxyz,qpm(1,:),qpm(1,:),qpm(2,:),11,qpmijk(1,:),xyz );
                       [ptsp3,tmpbrp3]=intersectfr( Nxyz,qpmstu(5,:),qpm(5,:),r,fn,ti,h,qpmijk(5,:),delta,arctheta,dirbr,xyz,btotal);
                       dstp=zeros(1,ptsp3);
                       for ij=1:ptsp3
                           tmpdir=tmpbrp3(ij,:)-r;
                           if ((dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir))<arctheta2)
                               dstp(ij)=norm(tmpdir);
                           elseif ((dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir))>=arctheta2)
                               dstp(ij)=0;
                               tmpbrp3(ij,:)=r;
                           end
                       end
                       lastd=repmat(d2,1,ptsp3);
                       [~,nnp3]=min(abs(dstp-lastd));
                       dnum=dnum+1;
                   end
                   if isempty(nnp3)==0
                       d2p3=norm(tmpbrp3(nnp3,:)-r);
                       crementp3=d2p3-d2p2;
                       if d2p2>deltam&&crementp3<eps||d2p2<deltap&&crementp3>-eps
                           brs=tmpbr;
                           di=d2;
                           bnn=nn;
                           continue;
                       end
                   end
                   
               end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                
                if d2p2>deltam&&pts>1&&ptsp2==1&&m>1
                    if m<extao-1
                        [ ~,~,qpm,~,qpmijk,~,qpmstu ] = geninit( Nxyz,qm(m+1,:),qm(1,:),qm(m+2,:),300,qmijk(m,:),xyz );
                        [ptsp3,tmpbrp3]=intersectfr( Nxyz,qpmstu(5,:),qpm(5,:),r,fn,ti,h,qpmijk(5,:),delta,arctheta,dirbr,xyz,btotal);
                        [flag,accbr] = accuracy( ptsp3,tmpbrp3,r,delta,dirbr,deltao,arctheta );
                        if flag==1
                            br=accbr;
                            jishu=m-1;
                            return;
                        end
                        if ptsp3~=0
                            dstp3=zeros(1,ptsp3);
                            for ij=1:ptsp3
                                tmpdir=tmpbrp3(ij,:)-r;
                                if ((dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir))<arctheta2)
                                    dstp3(ij)=norm(tmpdir);
                                elseif ((dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir))>=arctheta2)
                                    dstp3(ij)=0;
                                    tmpbrp3(ij,:)=r;
                                end
                            end
                            lastd3=repmat(delta,1,ptsp3);
                            [~,nnp3]=min(abs(dstp3-lastd3));
                            d2p3=norm(tmpbrp3(nnp3,:)-r);
                            crementp3=d2p3-d2p2;
                            if (d2>deltam&&crementp<eps)||((d2<deltap&&crementp>-eps)&&d2p2>deltam&&crementp3>-eps)
                                jishu=m;
                                brs=tmpbr;
                                return;
                            end
                            if m>1&&((crementp<eps&&d2<deltap&&crementp3>-eps)||(crementp3>-eps&&d2p2>deltam&&d2>deltam))                                
                                jishu=m-1;
                                arctheta3=larctheta3;
                                nn=bnn;
                                exbr=brs(nn,:);
                                bk=bbk;
                                return;
                            end
                        end
                    end
					[ ~,~,qpm,~,qpmijk,~,qpmstu ] = geninit( Nxyz,qm(m-1,:),qm(1,:),qm(m,:),10,qmijk(m,:),xyz );
                    [ptspb,tmpbrpb]=intersectfr( Nxyz,qpmstu(7,:),qpm(7,:),r,fn,ti,h,qpmijk(7,:),delta,arctheta,dirbr,xyz,btotal);
                    [flag,accbr] = accuracy( ptspb,tmpbrpb,r,delta,dirbr,deltao,arctheta );
                    if flag==1
                        br=accbr;
                        jishu=m-1;
                        return;
                    end
                    dstp=zeros(1,ptspb);
                    for ij=1:ptspb
                        tmpdir=tmpbrpb(ij,:)-r;
                        if ((dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir))<arctheta2)
                            dstp(ij)=norm(tmpdir);
                        elseif ((dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir))>=arctheta2)
                            dstp(ij)=0;
                            tmpbrpb(ij,:)=r;
                        end
                    end
                    lastd=repmat(d2,1,ptspb);
                    [~,nnpb]=min(abs(dstp-lastd));
                    dnum=0;
                    while (ptspb~=pts||nnpb~=nn)&&dnum<20
                        [ ~,~,qpm,~,qpmijk,~,qpmstu ] = geninit( Nxyz,qpm(9,:),qpm(1,:),qpm(10,:),11,qpmijk(1,:),xyz );
                        [ptspb,tmpbrpb]=intersectfr( Nxyz,qpmstu(5,:),qpm(5,:),r,fn,ti,h,qpmijk(5,:),delta,arctheta,dirbr,xyz,btotal);
                        dstp=zeros(1,ptspb);
                        for ij=1:ptspb
                            tmpdir=tmpbrpb(ij,:)-r;
                            if ((dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir))<arctheta2)
                                dstp(ij)=norm(tmpdir);
                            elseif ((dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir))>=arctheta2)
                                dstp(ij)=0;
                                tmpbrpb(ij,:)=r;
                            end
                        end
                        lastd=repmat(d2,1,ptspb);
                        [~,nnpb]=min(abs(dstp-lastd));
                        dnum=dnum+1;
                    end
                    if ptspb>0%%%%%%
                        d2pb=norm(tmpbrpb(nnpb,:)-r);
                        crementpb=d2-d2pb;
                        if any(dstp>0)
                            if (d2p2>deltam&&crementpb<eps&&d2<deltap)||(d2>deltam&&crementp>-eps&&crementpb>-eps)
                                jishu=m-1;
                                arctheta3=larctheta3;
                                nn=bnn;
                                bk=bbk;
                                exbr=brs(nn,:);
                                return;
                            end
                            jishu=m;
                            brs=tmpbr;
                            return;
                        end
                    end
                end
                                
                if m==1
                    [ ~,~,qpm,~,qpmijk,~,qpmstu ] = geninit( Nxyz,qm(m+1,:),qm(1,:),qm(m+2,:),11,qmijk(m,:),xyz );
                    [ptsp3,tmpbrp3]=intersectfr( Nxyz,qpmstu(5,:),qpm(5,:),r,fn,ti,h,qpmijk(5,:),delta,arctheta,dirbr,xyz,btotal);
                    if ptsp3~=0
                        [flag,accbr] = accuracy( ptsp3,tmpbrp3,r,delta,dirbr,deltao,arctheta );
                        if flag==1
                            br=accbr;
                            jishu=m-1;
                            return;
                        end
                        dstp3=zeros(1,ptsp3);
                        for ij=1:ptsp3
                            tmpdir=tmpbrp3(ij,:)-r;
                            if ((dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir))<arctheta2)
                                dstp3(ij)=norm(tmpdir);
                            elseif ((dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir))>=arctheta2)
                                dstp3(ij)=0;
                                tmpbrp3(ij,:)=r;
                            end
                        end
                        deltas3=repmat(delta,1,ptsp3);
                        [~,nnp3]=min(abs(dstp3-deltas3));
                        dnum=0;
                        while (ptsp3~=ptsp2||nnp3~=nnp2)&&dnum<20
                            [ ~,~,qpm,~,qpmijk,~,qpmstu ] = geninit( Nxyz,qpm(1,:),qpm(1,:),qpm(2,:),11,qpmijk(1,:),xyz );
                            [ptsp3,tmpbrp3]=intersectfr( Nxyz,qpmstu(5,:),qpm(5,:),r,fn,ti,h,qpmijk(5,:),delta,arctheta,dirbr,xyz,btotal);
                            dstp=zeros(1,ptsp3);
                            for ij=1:ptsp3
                                tmpdir=tmpbrp3(ij,:)-r;
                                if ((dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir))<arctheta2)
                                    dstp(ij)=norm(tmpdir);
                                elseif ((dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir))>=arctheta2)
                                    dstp(ij)=0;
                                    tmpbrp3(ij,:)=r;
                                end
                            end
                            lastd=repmat(d2,1,ptsp3);
                            [~,nnp3]=min(abs(dstp-lastd));
                            dnum=dnum+1;
                        end
                        d2p3=norm(tmpbrp3(nnp3,:)-r);
                        crementp3=d2p3-d2p2;
                        
                        if (pts==ptsp3&&nn==nnp3)==0
                            if (crementp2>-eps&&crementp3>-eps&&d2p2>deltam)&&(ptsp3==ptsp2&&nnp3==nnp2)
                                jishu=m;
                                brs=tmpbr;
                                return;
                            end
                        end
                    end
                end
                
                if m<extao
                    if d2p2<deltap&&crementp2>-eps||d2p2>deltam&&crementp2<eps
                        brs=tmpbr;
                        di=d2;
                        bnn=nn;
                        continue;
                    end
                end
                
                if cg==1&&((d2p2>di&&d2p2<deltap&&di<deltap&&nnp2<nnp)||(d2p2<di&&d2p2>deltam&&di>deltam&&nnp2>nnp))
                    brs=tmpbr;
                    di=d2;
                    bnn=nn;
                    continue;
                else
                    if crementp<eps&&d2>deltam
                        if (d2p2>d2||d2p2<deltap)&&cg==0
                            jishu=m;
                            brs=tmpbr;
                            return;
                        end
                    end
                    if ptsp==ptsp2&&((crementp<eps&&crementp2<eps&&d2p2<deltap&&d2p>deltam)||(crementp>-eps&&crementp2>-eps&&d2p2>deltam&&d2p<deltap))
                        brs=tmpbr;
                        di=d2;
                        bnn=nn;
                        continue;
                    end
                end
                
            end
            
            if ((crementp>-eps&&d2>deltam)||(crementp<eps&&d2<deltap))&&m>1
                jishu=m-1;
                arctheta3=larctheta3;
                nn=bnn;
                bk=bbk;
                exbr=brs(nn,:);
                return;
            end
        end
          
        crement=d2-di;  
        
        if m>1&&m<extao&&crement<eps&&pts==1
            [ptstr,tmpbrtr,~]=intersectfr( Nxyz,qmstu(m+1,:),qm(m+1,:),r,fn,ti,h,qmijk(m+1,:),delta,arctheta,dirbr,xyz,btotal);
            [flag,accbr] = accuracy( ptstr,tmpbrtr,r,delta,dirbr,deltao,arctheta );
            if flag==1
                br=accbr;
                jishu=m-1;
                return;
            end
            if ptstr==1
                 tmpdir=tmpbrtr-r;
                 if (dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir))>arctheta
                     brs=tmpbr;
                     di=d2;
                     bnn=nn;
                     continue;
                 end
            end
        end

        if crement>-eps&&di<deltap&&chgflg==0%if no intersection in admissable range which all |br-r|>delta 
            if d2>deltam
                jishu=m-1;
                arctheta3=larctheta3;
                nn=bnn;
                bk=bbk;
                if abs(md-delta)<abs(norm(brs(nn,:)-r)-delta)
                    exbr=mbr;
%                     jishu=mm;
%                     bk=mbk;
                else
                    if bbk==1
                        exbr=r;
                    else
                        exbr=brs(nn,:);
                    end
                end
                return;
            end      
        end
        if crement<eps&&di>deltam&&chgflg==0
            if d2<deltap
                jishu=m-1;
                arctheta3=larctheta3;
                nn=bnn;
                bk=bbk;
                if abs(md-delta)<abs(norm(brs(nn,:)-r)-delta)
                    exbr=mbr;
                else
                    if bbk==1
                        exbr=r;
                    else
                        exbr=brs(nn,:);
                    end
                end
                return;
            end
        end

        brs=tmpbr;
        
        di=d2;
        bnn=nn;
    else
        di=0;
        if bbk==0&&(m>ceil(extao/3)&&m>6)&&oflg==1&&dflg==1
            if di<deltam*0.2
                bk=1;
                continue;
            end
            nn=bnn;
            if abs(md-delta)<abs(norm(brs(nn,:)-r)-delta)
                exbr=mbr;
            else
                exbr=brs(nn,:);
            end
            jishu=m-1;
            bk=bbk;
            return;
        end
        if m>1
            bk=1;
        end
    end
    
end

jishu=m;
exbr=mbr;
brs=tmpbr;
if d2<deltam&&abs(md-delta)<abs(d2-delta)&&md>0
    exbr=mbr;
%     jishu=mm;
%     bk=mbk;
end


end

