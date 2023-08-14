function [ index] = anulus( id,numcr)
%
%   

if id>0&&id<=numcr
    index=id;
elseif id<=0
    index=id+numcr;
elseif id>numcr
    index=id-numcr;
end




end

