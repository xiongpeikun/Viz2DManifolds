function [points, tbr  ] = shootm( ptp,a,id,r,fn )
%SHOOTM Summary of this function goes here
% shooting method

eps=3e-14;
points=0;  %number of intersection points
tbr=zeros(id,3);

too=0;


%for intersections
if id==1&&ptp<eps
    points=1;
    tbr=a;
    tbr=tbr(1,:);
    return;
end

for i=2:id    
    if too==1
        too=0;
        continue;
    end
    if ptp(i-1)*ptp(i)<eps && ptp(i-1)*ptp(i)>=0
        if (ptp(i-1)<eps && ptp(i-1)>=0)||(ptp(i-1)>-eps && ptp(i-1)<=0)
            points=points+1;
            tbr(points,:)=a(i-1,:);
        end
        if (ptp(i)<eps && ptp(i)>=0)||(ptp(i)>-eps && ptp(i)<=0)
            points=points+1;
            tbr(points,:)=a(i,:);
            too=1;
        end
    elseif ptp(i-1)*ptp(i)<0
        points=points+1;
        pt1=a(i-1,:);
        pt2=a(i,:);
        d=(r-pt1)*fn'/norm(fn);
        p=pt1+d*(pt2-pt1)/norm(pt2-pt1);
        tbr(points,:)= p;
    end
   
end

if points>0
    tbr=tbr(1:points,:);
else
    tbr=[];
end

% if points==1
%     tbr=tbr(1,:);
% elseif points>1
%     tbr=tbr(1:points,:);
%     pnum=0;
%     inbr=zeros(points,3);
%     for i=2:points
%         if norm(tbr(i-1,:)-tbr(i,:))<eps
%            tbr(i,:)=tbr(i-1,:);
%            continue;
%         end
%         pnum=pnum+1;
%         inbr(pnum,:)=tbr(i-1,:);
%         if i==points
%             pnum=pnum+1;
%             inbr(pnum,:)=tbr(i,:);
%         end
%     end
%     inbr=inbr(1:pnum,:);
%     points=pnum;
%     tbr=inbr;
% else
%     tbr=[];
% end

end

