function v=volum(phase,P,T,Rg,amix,bmix)
e0=1+2^0.5;
s0=1-2^0.5;
c1=1;
c2=e0*bmix+s0*bmix-bmix-Rg*T./P;
c3=e0*s0*(bmix.^2)-(bmix.^2).*(e0+s0)-(e0+s0)*Rg*T*bmix./P+amix./P;
c4=-e0*s0*(bmix.^3)-e0*s0*Rg*T*(bmix.^2)./P-amix*bmix./P;
r1=roots([c1 c2 c3 c4]);
n1=0;
r2=[];
for i=1:3
    if isreal(r1(i))==1 & (r1(i)>0)
        r2(n1+1)=r1(i);
        n1=n1+1;
    end
end
r3=[];
r3=sort(r2,'descend');
if phase==1
    v=min(r3);
end
if phase==0
    v=max(r3);
end
end


