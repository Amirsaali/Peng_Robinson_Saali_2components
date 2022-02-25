function object=objective(par,np,P,T,Tc,Pc,Rg,mola,MW)
[Pcal Ycal]=pressure(par,np,P,T,Tc,Pc,Rg,mola,MW);
% [Leftt Rightt DelttaA]=consistency(par,np,P,T,Tc,Pc,Rg,mola,MW);
[Leftt Rightt DelttaA]=consistency(par,np,P,T,Tc,Pc,Rg,mola,MW);
object=0;

A=1; % weighting Factor for the objective Function
    for j=1:np
        object0(j)=abs((P(j)-Pcal(j)))*100/P(j);
object=A*sum(object+object0(j))%+((1-A)*sum(DelttaA'));
    end
par

[Rightt' Leftt'  DelttaA']
results=[T' P' mola' Pcal' object0']
maximum=max(object0)
mean=sum(object0)/np
sumDEL=sum(DelttaA);
arddelta=sum(DelttaA)/(np-1);
maxdel=max(DelttaA)
[maxdel arddelta]
end