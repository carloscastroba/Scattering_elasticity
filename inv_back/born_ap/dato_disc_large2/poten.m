function [q11,q12,q21,q22] = poten(z1,z2)

% potential Q=[q1 q2;q3 q4]. 
% It is written in column vector as [q1;q2;q3;q4]

ani = zeros(size(z1));
rr=sqrt(z1.^2+z2.^2);
ani(rr<0.8) = 1;
ani(rr<0.6) = 0;
ani(abs(z1)+abs(z2)<0.2) = 1.2;
%ani=max(0,0.9-rr.^2);

%rect=(real(z)>0).*(real(z)<1.2).*(imag(z)>0.1).*(imag(z)<0.4);
%rect=rect*0.4;

%tri=(real(z)+imag(z)<0).*(real(z)>-0.5).*(real(z)<-0.2).*(imag(z)>-0.4);
ani=ani*100;

q11=ani;
q12=ani;
q21=ani;
q22=ani;
% poten 2

%res = res + cos(pi*real(z));