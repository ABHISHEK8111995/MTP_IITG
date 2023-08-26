clc
clear all
x=0:0.01:10;
alpha=200;
f=1000;
zeta=0.7625;
g=76923; %mpa
den=7.850*10^-6; %kg/mm^3
d=30;l=1000;
vol=(pi/4)*d^2*l;
m=vol*den;
m1=5; R=40;
i=m1*R^2/2;
k=1;  %initial guess
di=inf;   
error=1; %in percentage
j=0;
r=0;
while(k>0 && di>2.0)   %di is decided by manufacturing possiblity of making hole in shaft
    q=sqrt(1+1/k^2);  %q=d0/di
    di=k*d;
    d0=q*di;
    j_s=pi*d^4/32;
    j_h=pi*(d0^4-di^4)/32;
    i_s=(i+m*d^2/24);
    i_h=(i+m*(d0^2+di^2)/24);
    j=j+1;  %for no.of iterations
    f_h=sqrt(g*j_h/(l*i_h));
    F=abs(sqrt(j_s/i_s)-sqrt(j_h/i_h))/sqrt(j_s/i_s);

    if (F*100<=error)     
        k=k-0.005;
        r=r+1;
        D0(r)=d0;
        Di(r)=di;
        freq(r)=f_h;            
    else
        k=k-0.005;
        
    end
    
end 

D=[D0;Di;freq]';
%          D1=[D0;Di;freq]';
%          D=flip(D1);
         p=length(D);
         q=[];
            for i1=1:p
                q(i1)=i1;
            end 
          a=[q',D];         
        [s,v] = listdlg('PromptString','Select the figures','Name','Figures Options',...
                'SelectionMode','single','ListSize',[250 300],...
                'ListString',num2str(a));
            
        if v==1
         j_h1=pi*(D0(s)^4-Di(s)^4)/32;   
         i_h1=(i+m*(D0(s)^2+Di(s)^2)/24);     
        else
            return
        end

 j_s1=pi*d^4/32;
 j_h1=pi*(D0(1)^4-Di(1)^4)/32;
 i_s1=(i+m*d^2/24);
 
 D0(s)
 Di(s)
 
f_s=sqrt(g*j_s1/(l*i_s1))
f_h=sqrt(g*j_h1/(l*i_h1))

% ws=f_s/(2*pi)
% wh=f_h/(2*pi)

% a=0:0.0001:0.01;
% sigma_s1=(3/8)*(alpha*f_s)*a.^2 + sqrt((f^2./(4*i_s1^2*f_s^2*a.^2))-zeta^2*f_s^2);
% sigma_s2=(3/8)*(alpha*f_s)*a.^2 - sqrt((f^2./(4*i_s1^2*f_s^2*a.^2))-zeta^2*f_s^2);

a=0.1:0.01:5;
sigma_s1=(3/8)*(alpha*a.^2/f_s) + sqrt((f^2./(4*f_s^2*a.^2))-zeta^2);
sigma_s2=(3/8)*(alpha*a.^2/f_s) - sqrt((f^2./(4*f_s^2*a.^2))-zeta^2);

sigma_h1=(3/8)*(alpha*a.^2/f_h) + sqrt((f^2./(4*f_h^2*a.^2))-zeta^2);
sigma_h2=(3/8)*(alpha*a.^2/f_h) - sqrt((f^2./(4*f_h^2*a.^2))-zeta^2);



% sigma_h1=(3/8)*(alpha*f_h)*a.^2+sqrt((f^2./(4*i_h1^2*f_h^2*a.^2))-zeta^2*f_h^2);
% sigma_h2=(3/8)*(alpha*f_h)*a.^2-sqrt((f^2./(4*i_h1^2*f_h^2*a.^2))-zeta^2*f_h^2);

figure(1)
plot(sigma_s1,a,'--g')
hold on
plot(sigma_s2,a,'--r')
hold on
plot(sigma_h1,a,'b',sigma_h2,a,'k')
grid minor

% figure(2)
% plot(sigma_h1,a,sigma_h2,a)
% grid minor














