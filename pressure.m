function [Pcal Ycal PHIL PHIG VL VG]=pressure(par,np,P,T,Tc,Pc,Rg,mola,MW)

for j=1:np
    %%%%% 2:gas    1:Ionic Liquid
    x(2)=mola(j)/(mola(j)+1000/MW);
    x(1)=1-x(2);
    t=T(j);
    p=P(j);
    [a b]=PRSEOS(par,t,Tc,Pc,Rg);
    l(1)=par(1);
    l(2)=par(2);
    tao=par(3);
    m=0;
    [amixl bmixl]=PRSmixEOS(l,tao,m,t,a,b,x);
    [dal dbl]=derivative(a,b,amixl,bmixl,t,l,tao,m,x);
    err=1;
    e0=1e-6;
    iter=0;
    ycal(2)=0.9;ycal(1)=1-ycal(2);
    while err>e0
        iter=iter+1;
        vl=volum(1,p,t,Rg,amixl,bmixl);
        phil=fugacity(vl,p,t,bmixl,amixl,Rg,dal,dbl);
        [amixv bmixv]=PRSmixEOS(l,tao,m,t,a,b,ycal);
        [dav dbv]=derivative(a,b,amixv,bmixv,t,l,tao,m,ycal);
        vg=volum(0,p,t,Rg,amixv,bmixv);
        phig=fugacity(vg,p,t,bmixv,amixv,Rg,dav,dbv);
        kval=phil./phig;
        s1=kval(1)*x(1)+kval(2)*x(2);
        err2=1;
        while err2>e0
            ycal(1)=kval(1)*x(1)/s1;
            ycal(2)=kval(2)*x(2)/s1;
            [amixv bmixv]=PRSmixEOS(l,tao,m,t,a,b,ycal);
            [dav dbv]=derivative(a,b,amixv,bmixv,t,l,tao,m,ycal);
            vg=volum(0,p,t,Rg,amixv,bmixv);
            phig=fugacity(vg,p,t,bmixv,amixv,Rg,dav,dbv);
            kval=phil./phig;
            s2=kval(1)*x(1)+kval(2)*x(2);
            err2=abs(s2-s1);
            s1=s2;
        end
        if abs(s1-1)<1e-6
            err=e0/10;
            Pcal(j)=p;
            Ycal(j,:)=ycal;
            PHIL(j,:)=phil;
            PHIG(j,:)=phig;
            VL(j,:)=vl;
            VG(j,:)=vg;
        else
            p=p*s1;
        end
    end
end
Pcal;
Ycal;
end

