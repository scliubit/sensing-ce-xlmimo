clear all
close all
%%

aperture_list = linspace(0.5,5,10);
implementations = 2000;
beta=4;
load(['../results/error_HOLO_' num2str(implementations) '_' num2str(beta) '.mat'])
RMSE_HOLO = sqrt(error_HOLO/implementations);
hold on
plot(aperture_list,RMSE_HOLO(1,:),'-o',LineWidth=2,MarkerSize=10,Color="#0072BD")
plot(aperture_list,RMSE_HOLO(2,:),'-s',LineWidth=2,MarkerSize=10,Color="#D95319")
plot(aperture_list,RMSE_HOLO(3,:),'-v',LineWidth=2,MarkerSize=10,Color="#EDB120")
plot(aperture_list,RMSE_HOLO(4,:),'-*',LineWidth=2,MarkerSize=10,Color="#7E2F8E")


clear error_HOLO
beta=2;
load(['../results/error_HOLO_' num2str(implementations) '_' num2str(beta) '.mat'])
RMSE_HOLO = sqrt(error_HOLO/implementations);
hold on
plot(aperture_list,RMSE_HOLO(1,:),':o',LineWidth=2,MarkerSize=10,Color="#0072BD")
plot(aperture_list,RMSE_HOLO(2,:),':s',LineWidth=2,MarkerSize=10,Color="#D95319")
plot(aperture_list,RMSE_HOLO(3,:),':v',LineWidth=2,MarkerSize=10,Color="#EDB120")
plot(aperture_list,RMSE_HOLO(4,:),':*',LineWidth=2,MarkerSize=10,Color="#7E2F8E")


clear error_HOLO
beta=1;
load(['../results/error_HOLO_' num2str(implementations) '_' num2str(beta) '.mat'])
RMSE_HOLO = sqrt(error_HOLO/implementations);
hold on
plot(aperture_list,RMSE_HOLO(1,:),'-.o',LineWidth=2,MarkerSize=10,Color="#0072BD")
plot(aperture_list,RMSE_HOLO(2,:),'-.s',LineWidth=2,MarkerSize=10,Color="#D95319")
plot(aperture_list,RMSE_HOLO(3,:),'-.v',LineWidth=2,MarkerSize=10,Color="#EDB120")
plot(aperture_list,RMSE_HOLO(4,:),'-.*',LineWidth=2,MarkerSize=10,Color="#7E2F8E")


axis([1 4.5 0 0.8])
legend({'$28\,{\rm GHz}$, $\eta=3$','$39\,{\rm GHz}$, $\eta=3$','$47\,{\rm GHz}$, $\eta=3$','$60\,{\rm GHz}$, $\eta=3$','$28\,{\rm GHz}$, $\eta=2$','$39\,{\rm GHz}$, $\eta=2$','$47\,{\rm GHz}$, $\eta=2$','$60\,{\rm GHz}$, $\eta=2$','$28\,{\rm GHz}$, $\eta=1$','$39\,{\rm GHz}$, $\eta=1$','$47\,{\rm GHz}$, $\eta=1$','$60\,{\rm GHz}$, $\eta=1$'},'Interpreter','Latex')
xlabel('XL-MIMO Aperture $\lambda N_{\rm BS}/2$ [m]','Interpreter','Latex');
ylabel('Localization Error RMSE','Interpreter','Latex');
title('(b) Localization RMSE','Interpreter','Latex')
set(gca,'fontsize',14);
set(gca,'fontname','Times New Roman')
box on
grid on