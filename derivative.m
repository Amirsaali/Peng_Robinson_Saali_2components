function [da db]=derivative(a,b,amix,bmix,T,l,tao,m,x)
    da(1)=2*a(1)*x(1)+(2*x(2)*(1-tao)*(sqrt(a(1)*a(2))))-amix;
    da(2)=2*a(2)*x(2)+(2*x(1)*(1-tao)*(sqrt(a(1)*a(2))))-amix;
    db(1)=b(1);
    db(2)=b(2);
end