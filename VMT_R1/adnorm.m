function [ nad ] = adnorm(i,numcr,adcr)
%ADNORM Summary of this function goes here
%   Detailed explanation goes here

nad=zeros(6,3);

il=anulus( i-1,numcr);
ir=anulus( i+1,numcr);

adl=adcr(il,:)-adcr(i,:);
adr=adcr(ir,:)-adcr(i,:);

nad(1,:)=cross(adl,adr);
nad(4,:)=adcr(i,:);

il=anulus( i-2,numcr);
ir=anulus( i-1,numcr);
adl=adcr(il,:)-adcr(ir,:);
adr=adcr(i,:)-adcr(ir,:);

nad(2,:)=cross(adl,adr);
nad(5,:)=adcr(ir,:);

il=anulus( i+1,numcr);
ir=anulus( i+2,numcr);
adl=adcr(i,:)-adcr(il,:);
adr=adcr(ir,:)-adcr(il,:);
nad(3,:)=cross(adl,adr);
nad(6,:)=adcr(il,:);



end

