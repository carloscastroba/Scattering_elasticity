function [q11,q12,q21,q22] = poten(z1,z2)

% potential Q=[q1 q2;q3 q4]. 
% It is written in column vector as [q1;q2;q3;q4]

rr=abs(z1)+abs(z2);
ani=max(0,0.9-rr)*1;

%rect=(real(z)>0).*(real(z)<1.2).*(imag(z)>0.1).*(imag(z)<0.4);
%rect=rect*0.4;

%tri=(real(z)+imag(z)<0).*(real(z)>-0.5).*(real(z)<-0.2).*(imag(z)>-0.4);

q11=ani;
q12=ani;
q21=ani;
q22=ani;
% poten 2

%res = res + cos(pi*real(z));