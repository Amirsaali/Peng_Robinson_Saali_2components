function [a b]=PRSEOS(par,T,Tc,Pc,Rg)
for i=1:2
    ac(i)=0.457235*((Rg*Tc(i))^2)/Pc(i);
    Tr(i)=T/Tc(i);
    E(i)=1/Tr(i)-Tr(i);
    if i==1
          betL(1)=1.1223;
          betL(2)=0.5966;
          b(i)=(0.077796*Rg*Tc(i)/Pc(i));
        alpha(i)=(1+((betL(1)*(1-Tr(i))))+(3*betL(2)*((1-Tr(i))^3)));
    end
    if i==2
          betg(1)=0.4699;
          betg(2)=0.1521;
       b(i)=(0.077796*Rg*Tc(i)/Pc(i));
        alpha(i)=(1+((betg(1)*(1-Tr(i))))+(3*betg(2)*((1-Tr(i))^3)));
    end
    a(i)=ac(i)*alpha(i);
end
end
