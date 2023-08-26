%composite skins
clear all
clc
format long
% graphite fiber data
% Ef=230e9;
% vf=0.3;
% Gf=22e9;
% pf=1800;
% % Epoxy resin matrix data
% Em=3.4e9;
% vm=0.3;
% Gm=1.308e9;
% pm=1200;
%  fiber data
% Ef=input('youngs modulus of fiber: ');
% vf=input('poissons ratio :');
% Gf=input('shear modulus : ');
% pf=input('density of fiber material: ');
% % matrix data
% Em=input('youngs modulus matrix : ');
% Gm=input('shear modulus : ');
% pm=input('density of matrix material: ');
% vm=input('poissons ratio :');
% Vf=input('volume fraction: ');
Vf=0.7;
Vm=1-Vf;
% pc=pf*Vf+pm*Vm;
% E1=Em*Vm+Ef*Vf;
% E2=(Ef*Em)/(Ef*Vm+Em*Vf);
% v12=vm*Vm+vf*Vf;
% G12=(Gm*Gf)/(Vf*Gm+Vm*Gf);
E1=141.2e9;
E2=9.720e9;
v12=0.28;
v21=E2*v12/E1;
G12=5.530e9;
n=input('number of plies:');
% h=input('thickness of upper and lower laminate:');
h=0.004;
to=0;
A=[0,0,0;0,0,0;0,0,0];
B=[0,0,0;0,0,0;0,0,0];
D=[0,0,0;0,0,0;0,0,0];
% ABBD matrix
i=1;
while i<=n
     theta=input('enter angle of ply:');
%     t=input('thickess of ply:');
    t=0.004;
    hk1=-(h/2)+to;
    to=to+t;
    hk2=-(h/2)+to;
    S=[(1/E1),-(v21)/E2,0;-(v12)/E1,1/E2,0;0,0,1/G12];
    Q=inv(S);
    T=[(cos(theta))^2,(sin(theta))^2,2*sin(theta)*cos(theta);(sin(theta))^2,(cos(theta))^2,-2*sin(theta)*cos(theta);-sin(theta)*cos(theta),sin(theta)*cos(theta),(cos(theta))^2-(sin(theta))^2];
    R=[1,0,0;0,1,0;0,0,2];
    Qred=(T)\Q*R*T/(R);
    A=A+Qred.*(hk2-hk1);
    B=B+Qred.*(((hk2)^2)/2-((hk1)^2)/2);
    D=D+Qred.*(((hk2)^3)/3-((hk1)^3)/3);
    i=i+1;
end
A;
B;
D;
ABBD=[A,B;B,D];
t1=h;
t3=h;
% t2=input('enter thickness of core:');
t2=0.02;
% lb=input('enter length of the beam:');
lbeam=1.0;
n=input('enter no of elements:');
% df=input('enter degree of freedom per node :');
df=4;
Bm=input('enter magnetic flux density:');
l=lbeam/n;
% b=input('enter breadth of the beam');
b=0.025;
no=n+1;
i1=(b*t1^3)/12;
i2=(b*t2^3)/12;
i3=(b*t3^3)/12;
rho1=1536;
rho2=970;
rho3=1536;
A1=b*t1;
A2=b*t2;
A3=b*t3;
J1=rho1*i1;
J2=rho2*i2;
J3=rho3*i3;
if Bm==0
    G=0.7037e6;
     ita=0.3;
else if Bm==0.1
        G=4.9302e6;
        ita=0.27917;
    else if Bm==0.2
            G=6.0465e6;
            ita=0.30208;
        else if Bm==0.3
                G=7.3023e6;
                ita=0.31667;
            else if Bm==0.4
                    G=8e6;
                    ita=0.3125;
                else if Bm==0.5
                        G=8.5116e6;
                        ita=0.30416;
                    else if Bm==0.6
                            G=8.6977e6;
                            ita=0.3;
                        else if Bm==0.7
                                G=8.6977e6;
                                ita=0.28958;
                            else 
                                G=8.698e6;
                                ita=0.288;
                            end
                        end
                    end
                end
            end
        end
    end
end
Gc=G*(1+ita*sqrt(-1));
Q1=[ABBD(1,1),ABBD(1,4);
    ABBD(1,4),ABBD(4,4)];
Q2 =[ABBD(1,3),ABBD(1,2),ABBD(1,5),ABBD(1,6);
     ABBD(1,6),ABBD(1,5),ABBD(4,5),ABBD(4,6)];
P1=[ABBD(2,1),ABBD(2,4);
    ABBD(3,1),ABBD(3,4);
    ABBD(2,4),ABBD(5,4);
    ABBD(3,4),ABBD(6,4)];
P2=[ABBD(2,3),ABBD(2,2),ABBD(2,5),ABBD(2,6);
    ABBD(3,3),ABBD(3,2),ABBD(3,5),ABBD(3,6);
    ABBD(2,6),ABBD(2,5),ABBD(5,5),ABBD(5,6);
    ABBD(3,6),ABBD(3,5),ABBD(6,5),ABBD(6,6)];
Q1;
Qt=(Q1-Q2*inv(P2)*P1);
syms x
u1=[(1-x/l);0;0;0;x/l;0;0;0];
u3=[0;(1-x/l);0;0;0;x/l;0;0];
w=[0;0;(1-3*(x/l)^2+2*(x/l)^3);x*(1-x/l)^2;0;0;(3*(x/l)^2-2*(x/l)^3);x*((x/l)^2-x/l)];
m=rho1*A1*(u1)*(u1.')+J1*diff(w,'x')*diff(w.','x')+(rho1*A1+rho2*A2+rho3*A3)*(w)*(w.')+rho3*A3*(u3)*(u3.')+J3*diff(w,'x')*diff(w.','x')+rho2*A2*((u1+u3)/2+(t1*diff(w,'x')-t3*diff(w,'x'))/4)*((u1.'+u3.')/2+(t1*diff(w.','x')-t3*diff(w.','x'))/4)+J2*((u1-u3)/t2+(t1*diff(w,'x')+t3*diff(w,'x'))/(2*t2))*((u1.'-u3.')/t2+(t1*diff(w.','x')+t3*diff(w.','x'))/(2*t2));
f1=b*(Qt(1,1)*diff(u1,'x')*diff(u1.','x')+Qt(2,1)*diff(u1,'x')*diff(w.','x',2)+Qt(1,2)*diff(u1,'x')*diff(w.','x',2)+Qt(2,2)*diff(w,'x',2)*diff(w.','x',2));
f3=b*(Qt(1,1)*diff(u3,'x')*diff(u3.','x')+Qt(2,1)*diff(u3,'x')*diff(w.','x',2)+Qt(1,2)*diff(u3)*diff(w.','x',2)+Qt(2,2)*diff(w,'x',2)*diff(w.','x',2));
f2=Gc*(A2*((u1-u3)/t2+((t1+t3+2*t2)*diff(w,'x'))/(2*t2))*((u1.'-u3.')/t2+((t1+t3+2*t2)*diff(w.'))/(2*t2)));
K1=int(f1,0,l);
K2=int(f2,0,l);
K3=int(f3,0,l);
Ke=K1+K2+K3;
Me=int(m,0,l);
i=1;
while i<=no*df
    j=1;
    while j<=no*df
        gstiff(i,j)=0;
        gmass(i,j)=0;
        j=j+1;
    end
    i=i+1;
end
% connectivity matrix
i=1;
while i<=n
    j=1;
    while j<=2
        C(i,j)=i+j-1;
        j=j+1;
    end
    i=i+1;
end
C;
% Global stiffness matrix & Global mass matrix
i=1;
while i<=n
    j=1;
    while j<=2
      p=C(i,j);
       k=1;
      while k<=2
        q=C(i,k);
         k1=1;
        while k1<=df
          k2=1;
          while k2<=df
             gstiff((p-1)*df+k1,(q-1)*df+k2)= gstiff((p-1)*df+k1,(q-1)*df+k2)+Ke((j-1)*df+k1,(k-1)*df+k2);
             gmass((p-1)*df+k1,(q-1)*df+k2)= gmass((p-1)*df+k1,(q-1)*df+k2)+Me((j-1)*df+k1,(k-1)*df+k2);
              k2=k2+1;
          end
        k1=k1+1;
        end 
      k=k+1;
     end
    j=j+1;
   end
  i=i+1;
end
bc=input('enter 1 for cantilever and 2 for simply supported 3 for Fixed-Fixed boundary conditons : ');
if bc==1
% Cantilever boundary conditions
gstiff(1:4,:)=[];
gstiff(:,1:4)=[];
gmass(1:4,:)=[];
gmass(:,1:4)=[];
% Simply supported boundary conditions
  else if bc==2
gstiff(no*df-1,:)=[];
gstiff(:,no*df-1)=[];
gstiff(no*2+3,:)=[];
gstiff(:,no*2+3)=[];
gstiff(3,:)=[];
gstiff(:,3)=[];
gmass(no*df-1,:)=[];
gmass(:,no*df-1)=[];
gmass(:,no*2+3)=[];
gmass(no*2+3,:)=[];
gmass(3,:)=[];
gmass(:,3)=[];
    else 
% Fixed-Fixed boundary conditions            
gstiff(no*df-3:no*df,:)=[];
gstiff(:,no*df-3:no*df)=[];
gmass(no*df-3:no*df,:)=[];
gmass(:,no*df-3:no*df)=[];
gstiff(1:4,:)=[];
gstiff(:,1:4)=[];
gmass(1:4,:)=[];
gmass(:,1:4)=[];
   end
end
fgmass=gmass;
fgstiff=gstiff;
Dfreq=inv(fgmass)*fgstiff; 
lamda=(eig(Dfreq));
W=(sqrt(lamda))
f=W/(2*pi);
f1=sort(f)
loss=imag(f1)./real(f1)
[mode,sqrW]=eig(Dfreq);
% mode shapes
Mode=real(mode);
Modes=zeros(no);
j=1;
while j<=no*df-4
    i=2;
    while i<=no
        Modes(i,j)=Mode(4*i-5,j);
        i=i+1;
    end
    j=j+1;
end
Modes;
lo=0:lb/n:lb;
plot(lo,Modes(:,no*df-4))
hold on
plot(lo,Modes(:,no*df-5))
hold on
plot(lo,Modes(:,no*df-6))
hold on
plot(lo,Modes(:,no*df-7))
hold on
plot(lo,Modes(:,no*df-8),'--')
hold on
plot(lo,Modes(:,no*df-9),'^')
legend('first','second','third','fourth','fifth')
sqrw1=real(sqrW);
sqrw2=imag(sqrW);
Wmat=sqrt(sqrW);
wmat=real(Wmat);
sqrwc=inv(wmat)*sqrw2;
% Finding time response of free vibration using Newmark beta Method
% c=sqrwc(no*df-5,no*df-5); 
% L=sqrw1(no*df-5,no*df-5); 
% N=1;
% % X0=input('enter initial displacements :');
% X0=0;
% % dX0=input('enter initial velocity :');
% dX0=0.5;
% % d2X0=input('enter initial acceleration :');
% d2X0=0;
% dt=input('enter time step :');
% alpha=input('enter alpha value :');
% beta=input('enter beta value :');
% a0=1/(beta*(dt^2));
% a1=alpha/(beta*(dt));
% a2=1/(beta*(dt));
% a3=1/(2*beta)-1;
% a4=alpha/beta-1;
% a5=(dt/2)*(alpha/beta-2);
% a6=dt*(1-beta);
% a7=beta*dt;
% Keff=L+a0*N+a1*c;
% I=input('enter number of iterations : ');
% Xtot=zeros(I);
% dXtot=zeros(I);
% d2Xtot=zeros(I);
% Xtot(1:I-1,:)=[];
% dXtot(1:I-1,:)=[];
% d2Xtot(1:I-1,:)=[];
% Xtot=X0;
% dXtot=dX0;
% d2Xtot=d2X0;
% i=2;
% while i<=I
%     Feff=N*(a0*X0+a2*dX0+a3*d2X0)+c*(a1*X0+a4*dX0+a5*d2X0);
%     Xt=inv(Keff)*Feff;
%     Xtot(:,i)=Xt;
%     d2Xt=a0*(Xt-X0)-a2*dX0-a3*d2X0;
%     dXt=a1*(Xt-X0)-a4*dX0-a5*d2X0;
%     dXtot(:,i)=dXt;
%     d2Xtot(:,i)=d2Xt;
%     d2X0=d2Xt;
%     dX0=dXt;
%     X0=Xt;
%     i=i+1;
% end
% t=linspace(0,dt*I,I);
% X1=Xtot(1,1:I);
% figure(2)
% grid
% hold on
% plot(t,X1)
% xlabel('time')
% ylabel('transverse displacement at the tip (m) ')
% title('time response plot')

% Frequency Response for Forced vibration 
F=input('enter magnitude of force:');
syms we
if bc==1
i=1;
while i<=no*df-4
     if i==no*df-5
      wred(i,1)=1;
        else
        wred(i,1)=0;    
    end
    i=i+1;
end
end
Mres=fgmass;
Kres=fgstiff;
wred
Z=-we^2*Mres+Kres./((2*pi)^2);   %to gt nf in terms of Hz div by 2pi sqr
Zinv=inv(Z);
Wres=abs(Zinv)*wred*F;
Wres(no*df-5)
if bc==3
    i=1;
while i<=no*df-8
     if i==7
      wred(i,1)=1;
        else
        wred(i,1)=0;    
    end
    i=i+1;
end
Mres=fgmass;
Kres=fgstiff;
wred
Z=-we^2*Mres+Kres./((2*pi)^2);
Zinv=inv(Z);
Wres=abs(Zinv)*wred*F;
Wres(7)
end
if bc==2
    i=1;
while i<=no*df-2
     if i==10
      wred(i,1)=1;
        else
        wred(i,1)=0;    
    end
    i=i+1;
end
Mres=fgmass;
Kres=fgstiff;
wred
Z=-we^2*Mres+Kres./((2*pi)^2);
Zinv=inv(Z);
Wres=abs(Zinv)*wred*F;
Wres(10)
end





