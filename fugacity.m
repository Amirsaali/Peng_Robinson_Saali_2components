function phi=fugacity(v,P,T,b,a,Rg,da,db)
e0=1+2^0.5;
s0=1-2^0.5;
for i=1:2
Z=P*v/(Rg*T);
B=b*P/(Rg*T);
I=(1/(s0-e0))*log((Z+(s0*B))/(Z+(e0*B)));
q=a/(b*Rg*T);
dq=q*(1+(da(i)/a)-(db(i)/b));
phi(i)=exp(((db(i)/b)*(Z-1))-log(Z-B)-(dq*I));
end
