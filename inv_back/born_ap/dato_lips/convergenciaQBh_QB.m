% Convergencia a la aproximación de born Q_B^h -> Q_B
load exp1_bs_M8.mat
NN=N;
x=x1r;
y=x2r;
QB8=real(QB_11);
%surf(x,y,QB)
load exp1_bs_M7.mat
x1rt=pagetranspose(x1r);
x2rt=pagetranspose(x2r);
QB7=griddedInterpolant(x1rt,x2rt,real(QB_11));
a(3) = 1/NN^2*sum(sum(abs(QB8-QB7(y,x)).^2));

load exp1_bs_M6.mat
x1rt=pagetranspose(x1r);
x2rt=pagetranspose(x2r);
QB6=griddedInterpolant(x1rt,x2rt,real(QB_11));
a(2)=1/NN^2*sum(sum(abs(QB8-QB6(y,x)).^2));

load exp1_bs_M5.mat
x1rt=pagetranspose(x1r);
x2rt=pagetranspose(x2r);
QB5=griddedInterpolant(x1rt,x2rt,real(QB_11));
a(1)=1/NN^2*sum(sum(abs(QB8-QB5(y,x)).^2));
sqrt(a)