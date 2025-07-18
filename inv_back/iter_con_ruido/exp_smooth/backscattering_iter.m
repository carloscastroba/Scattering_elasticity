%% Construcción de la aproximación de Born para el backscattering
clear all

%% Main parameters 
R=1;      % radius of a disc containing the support of V
Rz=2.1*R; % size of the physical domain [-Rz,Rz]x[-Rz,Rz]
M=6;      % number of points in the x-grid is N=2^M

%% Lame coef
lam=2;
mu=1;
rml=sqrt(2*mu+lam);
rm=sqrt(mu);
K=rml/rm;

%% Physical mesh
% Grid points: x
[x1,x2,N,h] = create_grid(M, Rz);

%% Fourier mesh
Mxi = M;
Nxi = 2^Mxi;  % same points in the physical and frequency spaces
Rxi = Nxi/(4*Rz); % radius of the xi-grid such that such xi-grid is 
                  % taken over the square [-Rxi,Rxi)^2. 
%  The relation Rxi=Nxi/(4*Rz) comes from the restriction that Nxi=N.
%  We can also compute Rz from Rxi with the same formula
[xi1,xi2,ndn,hxi] = create_grid(Mxi, Rxi);
C = 2^(Mxi-1)+1; % Index for locating the origin on the xi grid
xi1(C,C)=1.e-9; % Replace the origin with a near-zero value on the xi-grid
xi2(C,C)=1.e-9; % Replace the origin with a near-zero value on the xi-grid
absxi  = sqrt(xi1.^2+xi2.^2); % Módulo de xi

%% Matrix potential [ q11 q12 ; q21 q22 ]
[q11,q12,q21,q22] = poten(x1,x2); % this is a matrix here
[q11r,q12r,q21r,q22r] = poten(x1,x2); % this is a matrix here

%% Load scattering data 
load backs_data_M6.mat

% Grid points: x
[x1,x2,N,h] = create_grid(M, Rz);


% Aproximación de Born
 QBo_11=(1/(h^2))*ifftshift(ifftn(ifftshift(gorro_QB_11)));
 QBo_21=(1/(h^2))*ifftshift(ifftn(ifftshift(gorro_QB_21)));
 QBo_22=(1/(h^2))*ifftshift(ifftn(ifftshift(gorro_QB_22)));
 QBo_12=(1/(h^2))*ifftshift(ifftn(ifftshift(gorro_QB_12)));


%% Multiplicamos por funcion de soporte compacto
rr=sqrt(x1.^2+x2.^2);
cutoff=(rr<0.9);
QBo_11=QBo_11.*cutoff;
QBo_12=QBo_12.*cutoff;
QBo_21=QBo_21.*cutoff;
QBo_22=QBo_22.*cutoff;

 % pintamos la born approximation invirtiendo el campo lejano
figure(2)
clf
subplot(1,2,1)
mesh(x1,x2,real(QBo_11));
subplot(1,2,2)
mesh(x1,x2,imag(QBo_11))

figure(1)
clf
subplot(1,2,1)
mesh(x1,x2,real(q11r));
subplot(1,2,2)
mesh(x1,x2,imag(q11r))

%pause
 % Nueva iteracion
 QB_11=QBo_11;
 QB_21=QBo_21;
 QB_22=QBo_22;
 QB_12=QBo_12;

error_inf(1)=max(max(abs(real(q11r-QB_11))))+max(max(abs(real(q21r-QB_21))))+max(max(abs(real(q12r-QB_12))))...
    +max(max(abs(real(q22r-QB_22))));
error(1)=h^2*sum(sum(abs(real(q11r-QB_11))))+h^2*sum(sum(abs(real(q21r-QB_21))))+h^2*sum(sum(abs(real(q12r-QB_12))))...
    +h^2*sum(sum(abs(real(q22r-QB_22))));

for nn=1:19
q11=QB_11;
q12=QB_12;
q21=QB_21;
q22=QB_22;

gorro_QB_11=zeros(Nxi,Nxi);
gorro_QB_12=zeros(Nxi,Nxi);
gorro_QB_21=zeros(Nxi,Nxi);
gorro_QB_22=zeros(Nxi,Nxi);

for i=1:Nxi    
    parfor j=1:Nxi
        omega=pi*absxi(i,j);
        theta1=-xi1(i,j)/absxi(i,j);
        theta2=-xi2(i,j)/absxi(i,j);
        theta1o=-theta2;
        theta2o=theta1;
        % Definimos las 4 ondas incidentes (ui_pp,ui_ps) longitudinales
        % (ui_sp,ui_ss) transversales
        ui_pp_1=exp(1i*omega*(theta1*x1+theta2*x2))*theta1;
        ui_pp_2=exp(1i*omega*(theta1*x1+theta2*x2))*theta2;
        ui_ps_1=exp(1i*2*omega/(K+1)*(theta1*x1+theta2*x2))*theta1;
        ui_ps_2=exp(1i*2*omega/(K+1)*(theta1*x1+theta2*x2))*theta2;
        ui_ss_1=exp(1i*omega*(theta1*x1+theta2*x2))*theta1o;
        ui_ss_2=exp(1i*omega*(theta1*x1+theta2*x2))*theta2o;
        ui_sp_1=exp(1i*2*omega/(1/K+1)*(theta1*x1+theta2*x2))*theta1o;
        ui_sp_2=exp(1i*2*omega/(1/K+1)*(theta1*x1+theta2*x2))*theta2o;
        % Scattered solution (una para cada onda incidente)
        [w_pp_1,w_pp_2]=sol_LS_fun(lam,mu,rml*omega,x1,x2,Rz,R,M,ui_pp_1,ui_pp_2,q11,q12,q21,q22);
        [w_ps_1,w_ps_2]=sol_LS_fun(lam,mu,2*rml/(K+1)*omega,x1,x2,Rz,R,M,ui_ps_1,ui_ps_2,q11,q12,q21,q22);
        [w_ss_1,w_ss_2]=sol_LS_fun(lam,mu,rm*omega,x1,x2,Rz,R,M,ui_ss_1,ui_ss_2,q11,q12,q21,q22);
        [w_sp_1,w_sp_2]=sol_LS_fun(lam,mu,2*rm/(1/K+1)*omega,x1,x2,Rz,R,M,ui_sp_1,ui_sp_2,q11,q12,q21,q22);
        % Transformada de Fourier de Q theta evaluada en -2*omega theta
%         aux1=q11.*theta1+q12.*theta2;
%         gorro_Qtheta_1=h^2*sum(sum(exp(-1i*(-2)*omega*(theta1*x1+theta2*x2)).*aux1));
%         aux2=q21.*theta1+q22.*theta2;
%         gorro_Qtheta_2=h^2*sum(sum(exp(-1i*(-2)*omega*(theta1*x1+theta2*x2)).*aux2));
        % Amplitud pp multiplicada por 2*mu+lambda
        aux_pp_1=q11.*w_pp_1+q12.*w_pp_2;
        gorro_Qw_pp_1=h^2*sum(sum(exp(-1i*(-1)*omega*(theta1*x1+theta2*x2)).*aux_pp_1));
        aux_pp_2=q21.*w_pp_1+q22.*w_pp_2;
        gorro_Qw_pp_2=h^2*sum(sum(exp(-1i*(-1)*omega*(theta1*x1+theta2*x2)).*aux_pp_2));
        u_pp_infty_1=(gorro_Qw_pp_1.*theta1+gorro_Qw_pp_2.*theta2)*theta1;
        u_pp_infty_2=(gorro_Qw_pp_1.*theta1+gorro_Qw_pp_2.*theta2)*theta2;
        % Amplitud ps multiplicada por mu
        aux_ps_1=q11.*w_ps_1+q12.*w_ps_2;
        gorro_Qw_ps_1=h^2*sum(sum(exp(-1i*(-2*K/(K+1))*omega*(theta1*x1+theta2*x2)).*aux_ps_1));
        aux_ps_2=q21.*w_ps_1+q22.*w_ps_2;
        gorro_Qw_ps_2=h^2*sum(sum(exp(-1i*(-2*K/(K+1))*omega*(theta1*x1+theta2*x2)).*aux_ps_2));
        u_ps_infty_1=gorro_Qw_ps_1-(gorro_Qw_ps_1.*theta1+gorro_Qw_ps_2.*theta2)*theta1;
        u_ps_infty_2=gorro_Qw_ps_2-(gorro_Qw_ps_1.*theta1+gorro_Qw_ps_2.*theta2)*theta2;
        % Amplitud p
        u_p_infty_1=u_pp_infty_1+u_ps_infty_1;
        u_p_infty_2=u_pp_infty_2+u_ps_infty_2;
        % Transformada de Fourier de Q thetao evaluada en -2*omega theta
%         aux1o=q11.*theta1o+q12.*theta2o;
%         gorro_Qthetao_1=h^2*sum(sum(exp(-1i*(-2)*omega*(theta1*x1+theta2*x2)).*aux1o));
%         aux2o=q21.*theta1o+q22.*theta2o;
%         gorro_Qthetao_2=h^2*sum(sum(exp(-1i*(-2)*omega*(theta1*x1+theta2*x2)).*aux2o));
        % Amplitud ss multiplicada por mu
        aux_ss_1=q11.*w_ss_1+q12.*w_ss_2;
        gorro_Qw_ss_1=h^2*sum(sum(exp(-1i*(-1)*omega*(theta1*x1+theta2*x2)).*aux_ss_1));
        aux_ss_2=q21.*w_ss_1+q22.*w_ss_2;
        gorro_Qw_ss_2=h^2*sum(sum(exp(-1i*(-1)*omega*(theta1*x1+theta2*x2)).*aux_ss_2));
        u_ss_infty_1=gorro_Qw_ss_1-(gorro_Qw_ss_1.*theta1+gorro_Qw_ss_2.*theta2)*theta1;
        u_ss_infty_2=gorro_Qw_ss_2-(gorro_Qw_ss_1.*theta1+gorro_Qw_ss_2.*theta2)*theta2;
        % Amplitud sp multiplicada por 2*mu+lambda
        aux_sp_1=q11.*w_sp_1+q12.*w_sp_2;
        gorro_Qw_sp_1=h^2*sum(sum(exp(-1i*(-2/(K+1))*omega*(theta1*x1+theta2*x2)).*aux_sp_1));
        aux_sp_2=q21.*w_sp_1+q22.*w_sp_2;
        gorro_Qw_sp_2=h^2*sum(sum(exp(-1i*(-2/(K+1))*omega*(theta1*x1+theta2*x2)).*aux_sp_2));
        u_sp_infty_1=(gorro_Qw_sp_1.*theta1+gorro_Qw_sp_2.*theta2)*theta1;
        u_sp_infty_2=(gorro_Qw_sp_1.*theta1+gorro_Qw_sp_2.*theta2)*theta2;
        % Amplitud s
        u_s_infty_1=u_ss_infty_1+u_sp_infty_1;
        u_s_infty_2=u_ss_infty_2+u_sp_infty_2;
        % Transformada de Fourier de la aproximación de Born
        gorro_QB_11(i,j)=theta1*u_p_infty_1+theta1o*u_s_infty_1;
        gorro_QB_21(i,j)=theta1*u_p_infty_2+theta1o*u_s_infty_2;
        gorro_QB_22(i,j)=theta2*u_p_infty_2+theta2o*u_s_infty_2;
        gorro_QB_12(i,j)=theta2*u_p_infty_1+theta2o*u_s_infty_1;
    end
end
%save backs_data_M=5.mat

% Aproximación de Born
 QB_11=QBo_11-(1/(h^2))*ifftshift(ifftn(ifftshift(gorro_QB_11))).*cutoff;
 QB_21=QBo_21-(1/(h^2))*ifftshift(ifftn(ifftshift(gorro_QB_21))).*cutoff;
 QB_22=QBo_22-(1/(h^2))*ifftshift(ifftn(ifftshift(gorro_QB_22))).*cutoff;
 QB_12=QBo_12-(1/(h^2))*ifftshift(ifftn(ifftshift(gorro_QB_12))).*cutoff;
 
 

error_inf(nn+1)=max(max(abs(real(q11r-QB_11))))+max(max(abs(real(q21r-QB_21))))+max(max(abs(real(q12r-QB_12))))...
    +max(max(abs(real(q22r-QB_22))))
error(nn+1)=h^2*sum(sum(abs(real(q11r-QB_11))))+h^2*sum(sum(abs(real(q21r-QB_21))))+h^2*sum(sum(abs(real(q12r-QB_12))))...
    +h^2*sum(sum(abs(real(q22r-QB_22))))


%% Plot results
% pintamos la born approximation invirtiendo el campo lejano
figure(2)
clf
subplot(1,2,1)
mesh(x1,x2,real(QB_11));
subplot(1,2,2)
mesh(x1,x2,imag(QB_11))

figure(1)
clf
subplot(1,2,1)
mesh(x1,x2,real(q11r));
subplot(1,2,2)
mesh(x1,x2,imag(q11r))

%pause
end

save back_results_M6.mat


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Funciones
%%%%%%%%%%%%%%%%%%%%%%5

% Compute L-S solution
function [w1,w2]=sol_LS_fun(lam,mu,wvalue,x1,x2,Rz,R,M,ui1,ui2,q11,q12,q21,q22)
%% Main parameters 
% M  = number of points in the x-grid is N=2^M
% lam , mu = Lame coef
% wvalue = Energy
% (x1,x2) = Grid points
% (ui1,ui2) = incident wave 
% (w1,w2) = solution of the L-S equation

%% Compute Fourier coefficients of the Green function with kp and ks
kp=wvalue/sqrt(2*mu+lam);
ks=wvalue/sqrt(mu);

% Compute with the Green function
[Kj_11,Kj_12,Kj_21,Kj_22]=Fou_green_elas(mu,ks,kp,Rz,R,M);

%% Second term:  R_w (Q*u_i) 
% compute the fourier transform of the second hand term
ff1=fftshift(fftn(fftshift(q11.*ui1+q12.*ui2)));
ff2=fftshift(fftn(fftshift(q21.*ui1+q22.*ui2)));
% Multiply the matrix with the Green function by the FT of the second hand term
% and make the inverse fourier transform
f1 = fftshift(ifftn(fftshift(Kj_11.*ff1+Kj_12.*ff2)));
f2 = fftshift(ifftn(fftshift(Kj_21.*ff1+Kj_22.*ff2)));

% Reshape f as vertical vector of 2x2^(2M) components
f1 = f1(:);
f2 = f2(:);
f=[f1;f2];
%f=wvalue^2*[f1;f2];

%% Solve Lippmann-Schwinger equation u%sing gmres
%  The operator (I-R_w) is in the function GV_LS  
%w = GV_LS(f, Kj_11,Kj_12,Kj_21,Kj_22,q11,q12,q21,q22, M);
w = gmres('GV_LS', f, 10, 1e-7, 20, [], [], f, Kj_11,Kj_12,Kj_21,Kj_22,q11,q12,q21,q22, M);

% Reshape the two components [w1;w2] to original form
w1 = w(1:2^(2*M));
w2 = w(2^(2*M)+1:2^(2*M+1));

w1 = reshape(w1, 2^M, 2^M);
w2 = reshape(w2, 2^M, 2^M);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Generate mesh
function [z1,z2,N,h] = create_grid(m, s)
  % Creates matrices with meshpoints for a domain of length 2s with N points
  % at each dimension
    N = 2^m;
    h = 2*s/N;
    
    x = h*((-N/2):(N/2 - 1));
    [z1,z2] = meshgrid(x);
end

%%%%%%%%%%%%%%%%%%%%%%%%%

%% Fourier coef. Green function -(D^*+w^2I)^(-1)
function [Kj11,Kj12,Kj21,Kj22]=Fou_green_elas(mu,ks,kp,Rz,R,M)

[z1,z2,n,h] = create_grid(M, Rz);
rr=sqrt(z1.^2+z2.^2);
c = 2^(M-1)+1;
rr(c,c)=1.e-5; % avoid singularity at zero

tau=kp/ks;
F1=besselh(0,ks*rr)-(besselh(1,ks*rr)-tau*besselh(1,kp*rr))./(ks*rr);
F2=-besselh(0,ks*rr)+2*(besselh(1,ks*rr)-tau*besselh(1,kp*rr))./(ks*rr)...
    +tau^2*besselh(0,kp*rr);
G11=1i/(4*mu)*(F1+F2.*z1.^2./rr.^2);
G22=1i/(4*mu)*(F1+F2.*z2.^2./rr.^2);
G12=1i/(4*mu)*(F2.*z1.*z2./rr.^2);
G21=G12;

% Smooth truncation of Green's function near the boundary:
% We multiply the Green function by the smooth cutoff function eta(z) given by
% eta(z) = 1 if |z| < 2R,
% eta(z) = 1 - (|z|-2R)/(s-2R) if 2R < |z| < s,
% eta(z) = 0 if |z| > s.
s  = abs(min(min(z1)));
ep = s-2*R; % s > 2R
bigind       = rr>=s;
G11(bigind) = 0;
G22(bigind) = 0;
G12(bigind) = 0;
G21(bigind) = 0;
medind       = (rr<s) & (rr>2*R); % s > 2R
G11(medind) = G11(medind).*(1-(rr(medind)-2*R)/ep);
G22(medind) = G22(medind).*(1-(rr(medind)-2*R)/ep);
G12(medind) = G12(medind).*(1-(rr(medind)-2*R)/ep);
G21(medind) = G21(medind).*(1-(rr(medind)-2*R)/ep);

G11(c,c)=0;
G22(c,c)=0;
G12(c,c)=0;
G21(c,c)=0;

% Compute discrete Fourier transform of Green's function
Kj11 = -h^2*fftshift(fftn(fftshift(G11)));
Kj22 = -h^2*fftshift(fftn(fftshift(G22)));
Kj12 = -h^2*fftshift(fftn(fftshift(G12)));
Kj21 = -h^2*fftshift(fftn(fftshift(G21)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



