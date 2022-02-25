function [Leftt Rightt DelttaA]=consistency(par,np,P,T,Tc,Pc,Rg,mola,MW)
[Pcal Ycal PHIL PHIG VL VG]=pressure(par,np,P,T,Tc,Pc,Rg,mola,MW);
n1=length(T);
np=n1;
nt=n1;
for ii=1:np
    xx2=mola(ii)/(mola(ii)+1000/MW(1));
    xx1=1-xx2;
    tt=T(ii);
    pp=P(ii);
    vvl=VL(ii);
    PHIG2(ii)=PHIG(ii,2);
    PHIG1(ii)=PHIG(ii,1);
    PHIL2(ii)=PHIL(ii,2);
    PHIL1(ii)=PHIL(ii,1);
    yy1(ii)=Ycal(ii,1);
    yy2(ii)=Ycal(ii,2);
    ZZZ(ii)=pp*vvl/(Rg*tt);

A01(ii)=1/((ZZZ(ii)-1)*PHIL2(ii));
B01(ii)=(1-xx2)/(xx2*(ZZZ(ii)-1)*PHIL1(ii));
C01(ii)=1/(pp*xx2);
end
AA=trapz(PHIL2,A01);
BB=trapz(PHIL1,B01);
CC=trapz(P,C01);
Right0=AA+BB;
Left0=CC;
TCT0=100*abs((Right0-Left0)/(Left0));
TCT1=100*((Right0-Left0)/(Left0));
b=0;
while b <np-1
    b=b+1;
    ll2=PHIL(b:b+1,2);
    ll1=PHIL(b:b+1,1);
    hh=P(b:b+1)';
   
    for  c= b:b+1      
    xx(c,2)=mola(c)/(mola(c)+1000/MW(1));
    xx(c,1)=1-xx(c,2);
   K(c)=1/((ZZZ(c)-1)*(PHIL(c,2)));
   J(c)=(1-xx(c,2))/(xx(c,2)*(ZZZ(c)-1)*PHIL(c,1));
   H(c)=1/(P(c)*xx(c,2));
    end  

[K(b) K(b+1)];
AAA(b)=trapz(ll2,[K(b) K(b+1)]);
BBB(b)=trapz(ll1,[J(b) J(b+1)]);
CCC(b)=trapz(hh,[H(b) H(b+1)]);
    
if T(c)==T(c-1)
     Rightt(b)=AAA(b)+BBB(b);
Leftt(b)=CCC(b);
DelttaA(b)=100*abs((Rightt(b)-Leftt(b))/(Leftt(b)));
end
end