function [ flag,accbr ] = accuracy( numpts,interpts,r,delta,dirbr,deltao,arctheta )
% adapt new br within the condition
%   find br in the distance of delta


% deltao=1;
e=0.01;  %numerical tolerance
if delta>deltao
    deltam=delta-deltao+deltao*(1+e);
    deltap=delta-deltao+deltao*(1-e);
else
    deltam=delta*(1+e);
    deltap=delta*(1-e);
end

flag=0;
accbr=[];

for i=1:numpts
    arc=norm(interpts(i,:)-r);    %||br-r||
    tmpdir=interpts(i,:)-r;
    if (dirbr*tmpdir')/(norm(dirbr)*norm(tmpdir))<arctheta %if new point on the proper side
        if arc>=deltap&&arc<=deltam
            accbr=interpts(i,:);
            flag=1;
            return;
        end
    end
end




end

