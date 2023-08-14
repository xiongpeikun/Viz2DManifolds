function [ next,stu,ijk ] = whichcell( stu,ijk )
%WHICHCELL Summary of this function goes here
%   Detailed explanation goes here
next=0;
s=stu(1);t=stu(2);u=stu(3);
if s<0
    n=abs(floor(stu(1)/1));
    ijk(1)=ijk(1)-n;
    stu(1)=n+s;
    next=1;
end
if s>1
    n=abs(floor(stu(1)/1));
    ijk(1)=ijk(1)+n;
    stu(1)=s-n;
    next=1;
end
if t<0
    n=abs(floor(stu(2)/1));
    ijk(2)=ijk(2)-n;
    stu(2)=n+t;
    next=1;
end
if t>1
    n=abs(floor(stu(2)/1));
    ijk(2)=ijk(2)+n;
    stu(2)=t-n;
    next=1;
end
if u<0
    n=abs(floor(stu(3)/1));
    ijk(3)=ijk(3)-n;
    stu(3)=n+u;
    next=1;
end
if u>1
    n=abs(floor(stu(3)/1));
    ijk(3)=ijk(3)+n;
    stu(3)=u-n;
    next=1;
end
    

end

