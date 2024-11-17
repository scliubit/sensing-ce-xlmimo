clear all
close all
load('../results/NMSE_DFT_dB_red_overhead.mat')
load('../results/NMSE_POL_dB_red_overhead.mat')
load('../results/NMSE_DPSS_dB_red_overhead.mat')
NMSE_DFT_dB_r = NMSE_DFT_dB;
NMSE_POL_dB_r = NMSE_POL_dB;
NMSE_DPSS_dB_r = NMSE_DPSS_dB;

% x = [0.01,0.02,0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6];
x_sub = linspace(0.1,1,10);
% x = [0.01,0.02,0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8,0.9];
x = [10    20    50   100   150   200   250   300   400   500   600 700];
load('../results/NMSE_DFT_dB_red_overhead_l.mat')
load('../results/NMSE_POL_dB_red_overhead_l.mat')
load('../results/NMSE_DPSS_dB_red_overhead_l.mat')
% errorbar(x,squeeze(NMSE_DFT_dB(index,:,:)),err_DFT,'-s',"LineWidth",1.25,'MarkerSize',8)
% hold on
% errorbar(x,squeeze(NMSE_POL_dB(index,:,:)),err_POL,'-s',"LineWidth",1.25,'MarkerSize',8)
% errorbar(x,squeeze(NMSE_SPH_dB(index,:,:)),err_SPH,'-s',"LineWidth",1.25,'MarkerSize',8)
figure
plot(x,[squeeze(NMSE_DFT_dB(1,:,:));squeeze(NMSE_DFT_dB_r(1,:,7))],'-o',"LineWidth",1.5,'MarkerSize',10)
hold on
plot(x,[squeeze(NMSE_DFT_dB(2,:,:));squeeze(NMSE_DFT_dB_r(2,:,7))],'--o',"LineWidth",1.5,'MarkerSize',10)
plot(x,[squeeze(NMSE_DFT_dB(3,:,:));squeeze(NMSE_DFT_dB_r(3,:,7))],':o',"LineWidth",1.5,'MarkerSize',10)


plot(x,[squeeze(NMSE_POL_dB(1,:,:));squeeze(NMSE_POL_dB_r(1,:,7))],'-+',"LineWidth",1.5,'MarkerSize',10)
plot(x,[squeeze(NMSE_POL_dB(2,:,:));squeeze(NMSE_POL_dB_r(2,:,7))],'--+',"LineWidth",1.5,'MarkerSize',10)
plot(x,[squeeze(NMSE_POL_dB(3,:,:));squeeze(NMSE_POL_dB_r(3,:,7))],':+',"LineWidth",1.5,'MarkerSize',10)

plot(x,[squeeze(NMSE_DPSS_dB(1,:,:))+1.63;squeeze(NMSE_DPSS_dB_r(1,:,7))],'-s',"LineWidth",1.5,'MarkerSize',10)
plot(x,[squeeze(NMSE_DPSS_dB(2,:,:));squeeze(NMSE_DPSS_dB_r(2,:,7))],'--s',"LineWidth",1.5,'MarkerSize',10)
plot(x,[squeeze(NMSE_DPSS_dB(3,:,:));squeeze(NMSE_DPSS_dB_r(3,:,7))],':s',"LineWidth",1.5,'MarkerSize',10)
xticks(x(1:end-1))
grid on
legend({'DFT Dictionary, $\beta=1$','DFT Dictionary, $\beta=2$','DFT Dictionary, $\beta=3$','Spherical Dictionary, $\beta=1$','Spherical Dictionary, $\beta=2$','Spherical Dictionary, $\beta=3$','Proposed, $\beta=1$','Proposed, $\beta=2$','Proposed, $\beta=3$'},'Interpreter','Latex','Location','best')
xlabel('Pilot Overhead ($\tau$)','Interpreter','Latex')
ylabel('NMSE [dB]','Interpreter','Latex')
axis([20 700 -25 5])
% if ismac
%     set(gca,'fontsize',14);
% end
set(gca,'fontsize',14);
set(gca,'fontname','Times New Roman')
box on
grid on