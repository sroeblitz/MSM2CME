%parameters from deterministic tristability
p1 = 0.5;
p2 = 1;
a1 = 15;
a2 = 5; 
k1 = 1;
k2 = 1;
b1 = 0.05;
b2 = 0.05;
n1 = 6;
n2 = 6;
l1 = 1;
l2 = 1;
q1=5;
q2=1;
S1     = 4;
S2     = 4/12.9;    %scaling needed to correct for typo in previous file
%S1 = 2*a1*12.9/4;
%S1=15;
%q2=2;
%q1=1;
% q1=2.8169;
% q2=7.9352;
% S1=2.5648;
% k2=1.2086;

 
V=5e-12;    %Volume of macrophage in L
u=0.43*1e-11;  %typical concentration of X1,2 in mol/L

%--------------------------------------------------------
% Parameter transformation
nA = 6e+23; %Avrogado constant
uVnA = u*V*nA; %3e+2; conversion factor for a,b,k,p,S
a1 = a1*uVnA;
a2 = a2*uVnA;
b1 = b1*uVnA;
b2 = b2*uVnA;
k1 = k1*uVnA;
k2 = k2*uVnA;
p1 = p1*uVnA;
p2 = p2*uVnA;
%S1 = S1;
%S2 = S2; 