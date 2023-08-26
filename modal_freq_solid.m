function F1=modal_freq_solid(d,l,g,den,ip)
syms w
L11=l/10;
j11=pi*d^4/32;
% g=76923;
ma11=(den*j11*L11/6);

m=cell(1,10);
k=cell(1,10);

for i=1:10
    m{i}=ma11*[2 1;1 2];
    k{i}=(g*j11/L11)*[1 -1;-1 1];
end
    
  m{10}=m{10}+[0 0;0 ip];


%global Mass matrix
MASS_MATRIX1=zeros(11,11);
for i=1:10
    MASS_MATRIX1(i:i+1,i:i+1)=MASS_MATRIX1(i:i+1,i:i+1)+m{i};
end
  
%global stiffness matrix
STIFFNESS_MATRIX1=zeros(11,11);
for i=1:10
    STIFFNESS_MATRIX1(i:i+1,i:i+1)=STIFFNESS_MATRIX1(i:i+1,i:i+1)+k{i};
end

%boundary condition
MASS_MATRIX1=MASS_MATRIX1(2:11,2:11);
STIFFNESS_MATRIX1=STIFFNESS_MATRIX1(2:11,2:11);


MASS_MATRIX1=vpa(-w^2*MASS_MATRIX1);
A22=(MASS_MATRIX1+STIFFNESS_MATRIX1);


%solving for frequency
A22=det(A22);
frequency22=solve(A22,w);
freq1=double(frequency22);
F1=min(freq1(freq1>0))

end



