% Programa para dibujar
load exp1_bs_M7.mat
figure(1)
plot(x1r(N/2,:),q11(N/2,:),LineWidth=2,LineStyle="--")
hold on
plot(x1r(N/2,:),real(QB_11(N/2,:)),LineWidth=2)
axis([-1 1 -0.5 1.2])
legend('Q_{11}','Q_{B{11}}')
set(gca,'fontsize',18)
xlabel('x_1')
saveas(gcf,'exp_lip_20','epsc')
