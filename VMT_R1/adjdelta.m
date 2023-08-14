function [ deltas,tenddelta,con,convergence] = adjdelta( cr,adcr,deltas,taomaxs,deltamin )
%adjust delta and detect convergence


numcr=length(cr);
crt=cr';
convergence=0;          %boolean, 1 means convergence detected 
numconvergence=0;       %number of convergence detected
origin=[0,0,0];

con=zeros(1,numcr);
arcrb=zeros(1,numcr);
for i=1:numcr
    arcrb(i)=norm(adcr(i,:)-crt(i,:)); %record arclength of every leaf
    if arcrb(i)<deltamin*0.2
        numconvergence=numconvergence+1;
        con(i)=1;
    end
    if norm(origin-crt(i,:))<3
        con(i)=1;
    end
end

if numconvergence==numcr %if convergence detected on all leaves
    convergence=1;
end

cdeltas=deltas;
for i=1:numcr
    if taomaxs(i)==1
        cdeltas(i)=0;
    end
end


delta1=cdeltas(cdeltas>=deltamin);
[mind,~]=min(delta1);
if isempty(mind)
    [mind,~]=max(deltas);
end
if mind==0
    disp('all convergence');
    convergence=1;
end
tenddelta=mind;
for j=1:numcr
        
        if taomaxs(j)==1&&con(j)==0
            deltas(j)=mind;
        elseif con(j)==1
            deltas(j)=0;
        elseif isempty(mind)
             deltas(j)=deltamin;
        else
            deltas(j)=mind;
        end
    
end



end

