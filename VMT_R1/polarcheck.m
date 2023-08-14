function [ plusc,minusc ] = polarcheck( ilp,i,minumflg,numcr,tao )
%POLARCHECK Summary of this function goes here
%   Detailed explanation goes here

pplus=0;    
pminus=0;
if ilp>0&&ilp<=numcr
    if i>ilp%left
        for ii=ilp:i-1
            if minumflg(ii)==1
                pminus=pminus+1;
            end
            if minumflg(ii)==-1
                pplus=pplus+1;
            end
        end
    elseif i<ilp%right
        for ii=i+1:ilp
            if minumflg(ii)==1
                pminus=pminus+1;
            end
            if minumflg(ii)==-1
                pplus=pplus+1;
            end
        end
    end
elseif ilp<0% left
    il=ilp+numcr;
    for ii=il:numcr
        if minumflg(ii)==1
            pminus=pminus+1;
        end
        if minumflg(ii)==-1
            pplus=pplus+1;
        end
    end
    for ii=1:i-1
        if minumflg(ii)==1
            pminus=pminus+1;
        end
        if minumflg(ii)==-1
            pplus=pplus+1;
        end
    end
    
elseif ilp>numcr%right
    il=ilp-numcr;
    for ii=i+1:numcr
        if minumflg(ii)==1
            pminus=pminus+1;
        end
        if minumflg(ii)==-1
            pplus=pplus+1;
        end
    end
    for ii=1:il
        if minumflg(ii)==1
            pminus=pminus+1;
        end
        if minumflg(ii)==-1
            pplus=pplus+1;
        end
    end
end

plusc=pplus*(tao-1);
if pminus>0
    minusc=pminus*(tao-1)+1;
else
    minusc=0;
end


end

