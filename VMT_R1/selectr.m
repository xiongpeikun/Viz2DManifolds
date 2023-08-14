function [ r,rl,rr,ijkcube ] = selectr( i,numcr,cr,ijkindex )
%SELECTR 
%   

if i==1
    il=numcr;
    ir=i+1;
elseif i==numcr||i==0
    i=numcr;
    ir=1;
    il=i-1;
else
    il=i-1;
    ir=i+1;
end
r=cr(:,i)';                 % current r
rl=cr(:,il)';               % r_left
rr=cr(:,ir)';               % r_right
ijkcube=ijkindex(i,:);
% rlijk=ijkindex(il,:);
% rrijk=ijkindex(ir,:);
% ijkcube=[rijk;rlijk;rrijk];

end

