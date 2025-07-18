% Programa para dibujar
load exp2_bsM5.mat
figure(1)
plot(x1r(N/2,:),q11(N/2,:),LineWidth=2,LineStyle="--")
hold on
plot(x1r(N/2,:),real(QB_11(N/2,:)),LineWidth=2)
axis([-1 1 -0.5 2])
legend('q_{11}','q_{11}^B')
set(gca,'fontsize',18)
xlabel('x_1')
saveas(gcf,'prueba','epsc')
