clear all
close all
%%
implementations = 1000;

figure
% subplot(121)
hold on

load(['../results/error_MUSICc_' num2str(implementations) '.mat'])
load(['../results/error_OMP_' num2str(implementations) '.mat'])
load('../results/error_HOLO_2000.mat')
load(['../results/error_HOLOCc_' num2str(implementations) '.mat'])

beta = 1;
plot(sort(error_MUSIC(:,beta)),linspace(0,1,implementations),':','LineWidth',2,'Color',"#0072BD")
plot(sort(error_OMP(:,beta)),linspace(0,1,implementations),'--','LineWidth',2,'Color',"#0072BD")
plot(sort(error_HOLO(:,beta)),linspace(0,1,2000),'-','LineWidth',2,'Color',"#0072BD")
plot(sort(error_HOLO_C(:,beta)),linspace(0,1,implementations),'-.','LineWidth',2,'Color',"#0072BD")

beta = 2;
plot(sort(error_MUSIC(:,beta)),linspace(0,1,implementations),':','LineWidth',2,'Color',"#D95319")
plot(sort(error_OMP(:,beta)),linspace(0,1,implementations),'--','LineWidth',2,'Color',"#D95319")
plot(sort(error_HOLO(:,beta)),linspace(0,1,2000),'-','LineWidth',2,'Color',"#D95319")
plot(sort(error_HOLO_C(:,beta)),linspace(0,1,implementations),'-.','LineWidth',2,'Color',"#D95319")

beta = 3;
plot(sort(error_MUSIC(:,beta)),linspace(0,1,implementations),':','LineWidth',2,'Color',"#EDB120")
plot(sort(error_OMP(:,beta)),linspace(0,1,implementations),'--','LineWidth',2,'Color',"#EDB120")
plot(sort(error_HOLO(:,beta)),linspace(0,1,2000),'-','LineWidth',2,'Color',"#EDB120")
plot(sort(error_HOLO_C(:,beta)),linspace(0,1,implementations),'-.','LineWidth',2,'Color',"#EDB120")

legend({'MUSIC, $\eta=1$','BF, $\eta=1$','Proposed, $\eta=1$','STT, $\eta=1$','MUSIC, $\eta=2$','BF, $\eta=2$','Proposed, $\eta=2$','STT, $\eta=2$','MUSIC, $\eta=3$','BF, $\eta=3$','Proposed, $\eta=3$','STT, $\eta=3$'},'Interpreter','Latex',Location='best');

xlabel('Localization Error [m]','Interpreter','Latex');
ylabel('Localization Error CDF','Interpreter','Latex');
title('(a) Error CDF','Interpreter','Latex')
axis([0 0.7 0 1])
set(gca,'fontsize',14);
set(gca,'fontname','Times New Roman')
box on
grid on
