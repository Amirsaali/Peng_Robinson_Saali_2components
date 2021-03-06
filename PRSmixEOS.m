function [amix bmix]=PRSmixEOS(l,tao,m,T,a,b,x)
% if Classic==1
%     amix=a(1)*x(1)^2+a(2)*x(2)^2+2*(sqrt(a(1)*a(2)))*(1-l)*x(1)*x(2);
% %     amix=a(1)*x(1)^2+a(2)*x(2)^2+(sqrt(a(1)*a(2)))*(1-l1)*x(1)*x(2)+(sqrt(a(1)*a(2)))*(1-l2)*x(1)*x(2);
%     bmix=b(1)*x(1)+b(2)*x(2);
%   
% else
%     l(1)=par(5);
%     l(2)=par(6);
%     tao=par(7);
%     m=par(8);
%     l1=l(1);
%     l2=l(2);
%%%%%%%%%%%%%%%%%%%%% VDW 1 Mixing Rule %%%%%%%%%%%%%%%%%%%%
%     amix=a(1)*x(1)^2+a(2)*x(2)^2+(2*x(1)*x(2)*(1-l(1))*(sqrt(a(1)*a(2))));
%     bmix=b(1)*x(1)+b(2)*x(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% VDW 2 Mixing Rule %%%%%%%%%%%%%%%%%%%%
%     amix=a(1)*x(1)^2+a(2)*x(2)^2+(2*x(1)*x(2)*(1-l(1))*(sqrt(a(1)*a(2))));
%     bmix=b(1)*x(1)+b(2)*x(2)+(2*x(1)*x(2)*((b(1)+b(2))/2)*(1-tao)); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Shokouhi %%%%%%%%%%%%%%%%%%%% 
    amix=a(1)*x(1)^2+a(2)*x(2)^2+(2*x(1)*x(2)*(1-(x(1)*l(1))-(x(2)*l(2)))*(sqrt(a(1)*a(2))));
    bmix=b(1)*x(1)^2+b(2)*x(2)^2+(2*x(1)*x(2)*((b(1)+b(2))/2)*(1-tao));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end
end
