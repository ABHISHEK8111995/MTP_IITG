clc
clear all
format long g
% d=150;l=1000;
d=0.150;l=1;
% den=7.850*10^-6; %kg/mm^3
den=7829; %kg/mm^3
% g=76923; %mpa
e=206940*10^6;
u=0.288;
g=e/(2*(1+u))
vol=(pi/4)*d^2*l;
m=vol*den
m1=50; R=0.4
ip=m1*R^2/2
k=0.99;  %initial guess
% D=0;   
diff_freq=50; %in percentage  %variable
j=0;
r=0;
i_len=0.1*l;              %variable
F1=modal_freq_solid(d,l,g,den,ip);
while (i_len<0.5*l)   %variable    %di is decided by manufacturing possiblity of making hole in shaft
                                    %now it is 10 % of initial dia and initial length 
    while(k>0)                                      
        di=k*d;
        D=sqrt(d^2+(k^2*d^2*l)/i_len);      %d0=d
         F2=modal_freq_step(d,di,D,l,i_len,g,den,ip);
         F=abs(F1-F2)/F1;
        if (F*100<=diff_freq)
            k=k-0.5;            %variable
            r=r+1;
            D0(r)=D;
            Di(r)=di;
            x(r)=i_len;
            f(r)=F2;    %angular frequency
        else
            k=k-0.5;
        end
    end
    
    i_len=i_len+0.05*l;
    k=0.99;
%     D=0;
        
% dia_original=sqrt(((d^2-di^2)*(l-i_len)+(D^2-di^2)*i_len)/l); %this is mass conservation  
end 

D=[D0;Di;x;f];      % optimise_diameters



a=0.1:0.01:15;
alpha=200;
f=100000;
zeta=0.8;


sigma_s1=(3/8)*(alpha*a.^2/F1) + sqrt((f^2./(4*F1^2*a.^2))-zeta^2);
sigma_s2=(3/8)*(alpha*a.^2/F1) - sqrt((f^2./(4*F1^2*a.^2))-zeta^2);

sigma_h1=(3/8)*(alpha*a.^2/D(4)) + sqrt((f^2./(4*D(4)^2*a.^2))-zeta^2); %D(4) is first frquency of first iteration
sigma_h2=(3/8)*(alpha*a.^2/D(4)) - sqrt((f^2./(4*D(4)^2*a.^2))-zeta^2);

figure(1)
plot(sigma_s1,a,'--g')
hold on
plot(sigma_s2,a,'--r')
hold on
plot(sigma_h1,a,'b',sigma_h2,a,'k')
grid minor




