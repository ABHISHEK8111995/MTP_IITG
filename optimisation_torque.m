clc
clear all
x=0:0.1:10;
alpha=0.1;
g=79.3*10^3 %mpa
den=7850*10^-9; %kg/mm^3
d=150;l=1000;
vol=(pi/4)*d^2*l;
m=vol*den;
k=1; %initial guess 
i=225000; %mass=5 r=0.3
%d0=p*d this is limit
error=80; %in percentage
j=0;
z=0;
while(k>0)
%     q=sqrt(1+1/k^2);  %q=d0/di
%     di=k*d;
%     d0=d*k^4*((1+1/k^2)^2-1);
%     di=sqrt(d0^2/(1+1/k^2));
%     d0=q*di;
j=j+1;
q=[1 0 0 -1/k^3 -1];
sol=roots(q);
p1=sol(imag(sol)==0);
p2=p1(p1>0);
di=k*d;
d0=p2*di;
  j_s=pi*d^4/32;
  j_h=pi*(d0^4-di^4)/32;
  i_s=(i+m*d^2/24);
  i_h=(i+m*(d0^2+di^2)/24);
  z_s=pi*d^3/16;
  z_h=pi*(d0^4-di^4)/(16*d0);
  error_Freq=((j_s/i_s)-(j_h/i_h))/(j_s/i_s) 
  mass_s=pi*d^2/(4*l);
  mass_h=pi*(d0^2-di^2)/(4*l);
  error_mass=abs((mass_s-mass_h)/mass_s)*100;
   if (abs(error_Freq)*100<=error && error_mass<=error) 
%       disp(d0)
%       disp(di)
      area_solid=trapz(x,abs(g*j_s/l*ones(1,length(x))-g*j_s/l*(1+alpha*x.^2)));
      area_hollow=trapz(x,abs(g*j_h/l*ones(1,length(x))-g*j_h/l*(1+alpha*x.^2)));
      k=k-0.005;
      break
% %       if (area_solid>=area_hollow)
% %                  disp('THE OPTIMISED DIAMETERS ARE')
% %                  disp(di)
% %                  disp(d0)
% %                  disp(area_solid)
% %                  disp(area_hollow)
% %                   break
% %        end
%       
     
   else
       k=k-0.005;
   end
     
end 



k_s_l=(g*j_s/l).*x;
k_s_nl=(g*j_s/l*(1+alpha*x.^2)).*x;
k_hs_l=(g*j_h/l).*x;
k_hs_nl=(g*j_h/l*(1+alpha*x.^2)).*x;


figure(1)
plot(x,k_s_l,'b',x,k_s_nl,'--g')
hold on
plot(x,k_hs_l,'r',x,k_hs_nl,'--k')
xlabel('X')
ylabel('FORCE')
legend('solid-l','solid-nl','hollow-l','hollow-nl')
grid minor


figure(2)
plot(x,g*j_s/l*ones(1,length(x)),'b',x,g*j_s/l*(1+alpha*x.^2),'--g')
hold on
plot(x,g*j_h/l*(1+alpha*x.^2),'--k',x,(g*j_h/l)*ones(1,length(x)),'r')
xlabel('X')
ylabel('stifness')
legend('solid-l','solid-nl','hollow-nl','hollow-l')
grid minor








