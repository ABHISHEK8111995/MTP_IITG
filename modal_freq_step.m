function F2=modal_freq_step(d,di,D,l,i_len,g,den,ip)
% clc
% clear
% syms w
% d=150;
% di=15;
% D=160;
% l=1000;
% i_len=160;
% g=76923; %mpa
% den=7.850*10^-6; %kg/mm^3
% ip=2250000;
% l=50;

syms w
L1=(l-i_len)/5;
L2=(i_len)/5;

j1=pi*(d^4-di^4)/32;
j2=pi*(D^4-di^4)/32;
% g=76923;
ma1=(den*j1*L1/6);
ma2=(den*j2*L2/6);
% disp('............FINITE ELEMENT METHOD.........')

m=cell(1,10);
k=cell(1,10);

   for i=1:5
    m{i}=ma1*[2 1;1 2];
    k{i}=(g*j1/L1)*[1 -1;-1 1];
   end
   
    for i=6:10
      m{i}=ma2*[2 1;1 2];
      k{i}=(g*j2/L2)*[1 -1;-1 1];
    end
    
  m{10}=m{10}+[0 0;0 ip];

  
%global Mass matrix
MASS_MATRIX=zeros(11,11);
for i=1:10
    MASS_MATRIX(i:i+1,i:i+1)=MASS_MATRIX(i:i+1,i:i+1)+m{i};
end
  
%global stiffness matrix
STIFFNESS_MATRIX=zeros(11,11);
for i=1:10
    STIFFNESS_MATRIX(i:i+1,i:i+1)=STIFFNESS_MATRIX(i:i+1,i:i+1)+k{i};
end

%boundary condition
MASS_MATRIX=MASS_MATRIX(2:11,2:11);
STIFFNESS_MATRIX=STIFFNESS_MATRIX(2:11,2:11);


MASS_MATRIX=vpa(-w^2*MASS_MATRIX);
A2=(MASS_MATRIX+STIFFNESS_MATRIX);


%solving for frequency
A2=det(A2);
frequency2=solve(A2,w);
freq=double(frequency2);
F2=min(freq(freq>0));

end



