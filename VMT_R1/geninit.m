function [ fn,ql,qr,qlijk,qrijk,qlstu,qrstu ] = geninit( Nxyz,r,rl,rr,tao,ijkcube,xyz )
%GENINIT
% generating equally piecewise initial points on the front line

qlstu=zeros(tao,3);
qrstu=qlstu;
qlijk=qlstu;
qrijk=qlstu;
fn=rl-rr;
cdrrl=r-rl;                 % vector of r to r_left
cdrrr=r-rr;                 % vector of r to r_right


ql=repmat(r,tao,1);
qr=ql;


if tao==3
    ql(2,:)=(r+rl)/2;
    qr(2,:)=(r+rr)/2;
    ql(tao,:)=rl;
    qr(tao,:)=rr;
else
        
    morl=norm(cdrrl);
    morr=norm(cdrrr);
    stepmol=morl/(tao-1);
    stepmor=morr/(tao-1);
    for i=2:tao
        if i==tao
            ql(i,:)=rl;
            qr(i,:)=rr;
        else
            ql(i,:)=ql(i-1,:)-cdrrl/morl*stepmol;
            qr(i,:)=qr(i-1,:)-cdrrr/morr*stepmor;
        end
    end
end

rijk=ijkcube;
rcube=asigncube( xyz,rijk );


for m=1:tao                 %calculate initial points of qr(tao)
    %         m=1;
    [ nextl,stul,ijkl ] = insidecell( Nxyz,ql(m,:),rcube,rijk,xyz );
    [ nextr,stur,ijkr ] = insidecell( Nxyz,qr(m,:),rcube,rijk,xyz );
    if nextl==1
        qlstu(m,:)=stul;
        qlijk(m,:)=ijkl;
    elseif nextl==0
        qlijk(m,:)=rijk;
        qlstu(m,:)=xyztostu( Nxyz,ql(m,1),ql(m,2),ql(m,3),rijk,xyz );
    end
    if nextr==1
        qrstu(m,:)=stur;
        qrijk(m,:)=ijkr;
    elseif nextr==0
        qrijk(m,:)=rijk;
        qrstu(m,:)=xyztostu( Nxyz,qr(m,1),qr(m,2),qr(m,3),rijk,xyz );
    end
end
 

    



end

