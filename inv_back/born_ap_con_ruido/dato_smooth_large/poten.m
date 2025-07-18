function [q11,q12,q21,q22] = poten(z1,z2)

% potential Q=[q1 q2;q3 q4]. 
% It is written in column vector as [q1;q2;q3;q4]

smooth=max(0,2*exp(-5*((z1-0.5).^2+z2.^2)) +1.5*exp(-4*((z1+0.5).^2+(z2-0.4).^2))...
    +3*exp(-7*((z1+0.3).^2+(z2+0.3).^2))-0.8);
smooth=smooth.*((z1.^2+z2.^2)<1);

smooth=smooth*10;

q11=smooth;
q12=smooth;
q21=smooth;
q22=smooth;
% poten 2

%res = res + cos(pi*real(z));