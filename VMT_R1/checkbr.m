function [skipm] = checkbr(cr,br,numcr,deltafmin)
%CHECKBR Summary of this function goes here
%   Function of checking Mesh Concavities and Entanglements

mcr=[cr(:,end-1:end),cr,cr(:,1:2)];
mbr=[br(end-1:end,:);br;br(1:2,:)];
skipm=zeros(1,numcr+4);

for i=3:numcr+2
    mc1=(mcr(:,i-1)-mcr(:,i))/norm(mcr(:,i-1)-mcr(:,i));
    mc2=(mcr(:,i+1)-mcr(:,i))/norm(mcr(:,i+1)-mcr(:,i));
    n1=cross(mc1,mc2);
    dirc=dot(mc1,mc2);
    theta=dot(mbr(i,:)-mcr(:,i)',n1)/(norm(n1)*norm(mbr(i,:)-mcr(:,i)'));
    
    mb2=(mbr(i,:)-mbr(i-1,:))/norm(mbr(i,:)-mbr(i-1,:));
    
    mb4=(mbr(i,:)-mbr(i+1,:))/norm(mbr(i,:)-mbr(i+1,:));
    
    n2=cross(mb2,mb4);
    nm=dot(n1,n2);
    dir=mb2*mb4';

    if dirc>-0.1&&(abs(theta)<0.1||(isnan(theta)&&dir>dirc))
        if nm<0
            d1=norm(mbr(i-1,:)-mcr(:,i-1)');
            d2=norm(mbr(i+1,:)-mcr(:,i+1)');
            dc1=dot(mcr(:,i)-mcr(:,i+1),mbr(i-1,:)'-mcr(:,i+1))/norm(mcr(:,i)-mcr(:,i+1));
            dc2=dot(mcr(:,i)-mcr(:,i-1),mbr(i+1,:)'-mcr(:,i-1))/norm(mcr(:,i)-mcr(:,i-1));
            d1o=mcr(:,i)+(mcr(:,i+1)-mcr(:,i))/norm(mcr(:,i)-mcr(:,i+1))*dc1;
            d2o=mcr(:,i)+(mcr(:,i-1)-mcr(:,i))/norm(mcr(:,i)-mcr(:,i-1))*dc2;
            dis1=norm(d1o'-mbr(i-1,:));
            dis2=norm(d2o'-mbr(i+1,:));
            
            if d1>dis1||d2>dis2||isnan(theta)
                skipm(i)=1;
                skipm(i-1)=1;
                skipm(i+1)=1;
                if norm(mbr(i-2,:)-mbr(i+2,:))<deltafmin
                    skipm(i+2)=1;
                end
            end
        end
        if dir>0.5
            dbc=norm(mbr(i,:)-mcr(:,i)');
            dd1=dot(mcr(:,i)-mcr(:,i+1),mbr(i,:)'-mcr(:,i+1))/norm(mcr(:,i)-mcr(:,i+1));
            dd2=dot(mcr(:,i)-mcr(:,i-1),mbr(i,:)'-mcr(:,i-1))/norm(mcr(:,i)-mcr(:,i-1));
            d1o1=mcr(:,i)+(mcr(:,i+1)-mcr(:,i))/norm(mcr(:,i)-mcr(:,i+1))*dd1;
            d2o2=mcr(:,i)+(mcr(:,i-1)-mcr(:,i))/norm(mcr(:,i)-mcr(:,i-1))*dd2;
            diso1=norm(d1o1'-mbr(i,:));
            diso2=norm(d2o2'-mbr(i,:));
            
            nombdp1=(mbr(i,:)-mcr(:,i)')*((mcr(:,i+1)-mcr(:,i))/norm(mcr(:,i+1)-mcr(:,i)));
            bdp1=nombdp1*((mcr(:,i+1)-mcr(:,i))/norm(mcr(:,i+1)-mcr(:,i)));
            dp1=mcr(:,i)+bdp1;
            nombdm1=(mbr(i,:)-mcr(:,i)')*((mcr(:,i-1)-mcr(:,i))/norm(mcr(:,i-1)-mcr(:,i)));
            bdm1=nombdm1*((mcr(:,i-1)-mcr(:,i))/norm(mcr(:,i-1)-mcr(:,i)));
            dm1=mcr(:,i)+bdm1;
            
            touyingy=norm(mbr(i,:)-dp1');
            touyingz=norm(mbr(i,:)-dm1');
            
            if dbc>diso1||dbc>diso2||touyingy<dbc*0.95||touyingz<dbc*0.95
                skipm(i)=1;
            end
            if norm(mbr(i-1,:)-mbr(i+1,:))<deltafmin
                skipm(i-1)=1;
            end
        end
        
    end

    
end
for i=1:3
    if skipm(i)==1
        skipm(end-4+i)=skipm(i);
    end
end

skipm=skipm(3:end-2);


end

