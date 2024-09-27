function [ m, bk ] = polarities( Nxyz,initao,dflg2,interjishu,qstu,q,qijk,r,fn,ti,h,delta,deltam,arctheta,dirbr,xyz,btotal )
%POLARITIES Summary of this function goes here
%   Detailed explanation goes here
% ,minum [numcr,2]
arctheta2=0.99;
bk=1;
taolen=length(q(:,1));
dflg=dflg2(1);
tao=dflg2(2);
num=0;
% m=1+(interjishu-1)*(tao-1)+(initao-1);
m=interjishu-2;
if taolen<m&&taolen<41
    m=0;
    bk=0;
    return;
end
range=40+initao;

while m>0&&num<range%&&taolen>=m
    num=num+1;
    [pts,tbr,oflg]=intersectfr( Nxyz,qstu(m,:),q(m,:),r,fn,ti,h,qijk(m,:),delta,arctheta,dirbr,xyz,btotal);
    if pts==0
        if oflg==1&&dflg==1%&&bbk==1
            break;
        else
            m=m-(tao-1);
        end
    else
        dst=zeros(1,pts);
        for ij=1:pts
            tmpdir=tbr(ij,:)-r;
            if ((dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir))<arctheta2)
                dst(ij)=norm(tmpdir);
            elseif ((dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir))>=arctheta2)&&norm(tmpdir)>deltam
                dst(ij)=0;
                tbr(ij,:)=r;
            end
        end
        deltas=repmat(delta,1,pts);
        [~,nn]=min(abs(dst-deltas));

        tmpdir=tbr(nn,:)-r;
        disbr=norm(tmpdir);
        tmparc=(dirbr*tmpdir')/(norm(dirbr)*disbr);
        
        if m<initao
            m=0;
            bk=0;
            return;
        end
        
        if (tmparc<-0.9&&disbr>0)
            mi=m-(tao-1);
            [ptsi,tbri,~]=intersectfr( Nxyz,qstu(mi,:),q(mi,:),r,fn,ti,h,qijk(mi,:),delta,arctheta,dirbr,xyz,btotal);
            if ptsi>0
                dsti=zeros(1,ptsi);
                for ij=1:ptsi
                    tmpdir=tbri(ij,:)-r;
                    if ((dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir))<arctheta2)
                        dsti(ij)=norm(tmpdir);
                    elseif ((dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir))>=arctheta2)&&norm(tmpdir)>deltam
                        dsti(ij)=0;
                        tbri(ij,:)=r;
                    end
                end
                deltas=repmat(delta,1,ptsi);
                [~,nn]=min(abs(dsti-deltas));
                
                tmpdiri=tbri(nn,:)-r;
                disbri=norm(tmpdiri);
                tmparci=(dirbr*tmpdiri')/(norm(dirbr)*disbri);
                if isnan(tmparci)||tmparci>arctheta
                    m=mi-2;
                    return;
                end
                if disbri>deltam&&tmparci<arctheta
                    m=mi-(tao-1);
                    continue;
                end
            end
        end
        
        if (disbr>deltam&&tmparc>arctheta&&tmparc<0.8)
            mi=m-(tao-1);
            [ptsi,tbri,~]=intersectfr( Nxyz,qstu(mi,:),q(mi,:),r,fn,ti,h,qijk(mi,:),delta,arctheta,dirbr,xyz,btotal);
            if ptsi>0
                dsti=zeros(1,ptsi);
                for ij=1:ptsi
                    tmpdir=tbri(ij,:)-r;
                    if ((dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir))<arctheta2)
                        dsti(ij)=norm(tmpdir);
                    elseif ((dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir))>=arctheta2)&&norm(tmpdir)>deltam
                        dsti(ij)=0;
                        tbri(ij,:)=r;
                    end
                end
                deltas=repmat(delta,1,ptsi);
                [~,nn]=min(abs(dsti-deltas));
                
                tmpdiri=tbri(nn,:)-r;
                disbri=norm(tmpdiri);
                tmparci=(dirbr*tmpdiri')/(norm(dirbr)*disbri);
                
                if disbri>deltam&&tmparci<arctheta
                    m=mi-(tao-1);
                    continue;
                end
            end
        end
        
        if disbr>deltam&&tmparc<arctheta
            m=m-(tao-1);
        else
            m=m-2;
            if disbr<deltam*0.01&&tmparc<arctheta
                m=m-1;
            end
            break;
        end
        if m<initao
            m=0;
            bk=0;
            return;
        end
        
    end

end
  

if num==range||m<0%||taolen<m
    m=0;
    bk=0;
end


end

