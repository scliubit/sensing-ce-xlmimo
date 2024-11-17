clear all
close all
load('../results/NMSE_DFT_dB_135.mat')
load('../results/NMSE_POL_dB_135.mat')
load('../results/NMSE_DPSS_dB_135.mat')
NMSE_DFT_dB_135 = NMSE_DFT_dB;
NMSE_POL_dB_135 = NMSE_POL_dB;
NMSE_DPSS_dB_135 = NMSE_DPSS_dB;

OMP_iternum = 45;
x = 1:OMP_iternum;
x_sub = 1:2:OMP_iternum;
% errorbar(x,squeeze(NMSE_DFT_dB(index,:,:)),err_DFT,'-s',"LineWidth",1.25,'MarkerSize',8)
% hold on
% errorbar(x,squeeze(NMSE_POL_dB(index,:,:)),err_POL,'-s',"LineWidth",1.25,'MarkerSize',8)
% errorbar(x,squeeze(NMSE_SPH_dB(index,:,:)),err_SPH,'-s',"LineWidth",1.25,'MarkerSize',8)
figure
plot(x_sub,squeeze(NMSE_DFT_dB_135(1,:,x_sub)),'-o',"LineWidth",1.5,'MarkerSize',10)
hold on
plot(x_sub,squeeze(NMSE_DFT_dB_135(2,:,x_sub)),'--o',"LineWidth",1.5,'MarkerSize',10)
plot(x_sub,squeeze(NMSE_DFT_dB_135(3,:,x_sub)),':o',"LineWidth",1.5,'MarkerSize',10)


plot(x_sub,squeeze(NMSE_POL_dB_135(1,:,x_sub)),'-+',"LineWidth",1.5,'MarkerSize',10)
plot(x_sub,squeeze(NMSE_POL_dB_135(2,:,x_sub)),'--+',"LineWidth",1.5,'MarkerSize',10)
plot(x_sub,squeeze(NMSE_POL_dB_135(3,:,x_sub)),':+',"LineWidth",1.5,'MarkerSize',10)

plot(x_sub,smooth(squeeze(NMSE_DPSS_dB_135(1,:,x_sub)),7),'-s',"LineWidth",1.5,'MarkerSize',10)
plot(x_sub,smooth(squeeze(NMSE_DPSS_dB_135(2,:,x_sub)),7),'--s',"LineWidth",1.5,'MarkerSize',10)
plot(x_sub,smooth(squeeze(NMSE_DPSS_dB_135(3,:,x_sub)),7),':s',"LineWidth",1.5,'MarkerSize',10)

grid on
legend({'DFT Dictionary, $\mu=0.1$','DFT Dictionary, $\mu=0.3$','DFT Dictionary, $\mu=0.5$','Spherical Dictionary, $\mu=0.1$','Spherical Dictionary, $\mu=0.3$','Spherical Dictionary, $\mu=0.5$','Proposed, $\mu=0.1$','Proposed, $\mu=0.3$','Proposed, $\mu=0.5$'},'Interpreter','Latex','Location','best')
xlabel('Number of OMP Iterations ($\it I$)','Interpreter','Latex')
ylabel('NMSE[dB]','Interpreter','Latex')
axis([1 OMP_iternum-1 -25 0])
% if ismac
%     set(gca,'fontsize',14);
% end
set(gca,'fontsize',14);
set(gca,'fontname','Times New Roman')
box on
grid on