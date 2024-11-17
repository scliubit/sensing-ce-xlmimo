clear all
close all
clc

%% usually takes 1 hour

addpath('../utils/')
fc = 28e9;
c = physconst('Lightspeed');
lambda = c/fc;
N = 256;
N_UE = 4;
d = lambda/2;
aperture = (N-1)*d; % ~1m
k = 2*pi/lambda;
bandwidth = 9e6;
L=6;
Rician_factor = 1;%~10dB

NT_coord.x = linspace(-aperture/2,aperture/2,N);
NT_coord.y = zeros(1,N);
% H_downlink = zeros(1,N);
% UE location range
UE_range_x = 3;
UE_range_y_min = 3;
UE_range_y_max = 20;


beta_list = [1,2,3];
%% SIM start
epsilon = 0;
% activeRatio=[0.1 0.2 0.25 0.3 0.5 0.75 1];
sampleRatio=linspace(0.1,1,10);%[0.4]
sampleRatio = [0.01,0.02,0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6];
sampleNumList=round(sampleRatio*N*N_UE);
sampleNumList = [10    20    50   100   150   200   250   300   400   500   600];
OMP_iternum = 30;
implementations = 1000; % no larger than the reference error size
[aR_DPSS,~] = dpss(N_UE,k*aperture/4/pi/UE_range_y_max/N*N_UE,N_UE,"trace");
NMSE_DFT = zeros(numel(beta_list),implementations,numel(sampleNumList));
NMSE_POLAR = zeros(numel(beta_list),implementations,numel(sampleNumList));
NMSE_DPSS = zeros(numel(beta_list),implementations,numel(sampleNumList));
tic

for beta_i = 1:numel(beta_list)
    beta = beta_list(beta_i);
    load(['../results/error_HOLO_1000_' num2str(beta) '.mat'])

    [Psi_DFT,Psi_POL,Psi_DPSS,index2coor] = dict_design_red(N_UE, N,aperture,k,UE_range_y_min,UE_range_y_max,UE_range_x,beta);

    for imp=1:implementations
        UE.x=2*(rand()-0.5)*UE_range_x;
        UE.y=UE_range_y_min+(UE_range_y_max-UE_range_y_min)*rand();
        aperture_UE = (N_UE-1)*d;
        NR_coord.x = UE.x+linspace(-aperture_UE/2,aperture_UE/2,N_UE);
        NR_coord.y = UE.y+zeros(1,N_UE);
        if L>=2
            theta = linspace(-pi,pi,L-1);
            R = 15;
            scatter_coord.x = R*cos(theta);
            scatter_coord.y = R*sin(theta);
        else
            scatter_coord.x = [];
            scatter_coord.y = [];
        end
        %         H_Downlink = GenChannelDL(N_UE,N,L,NR_coord,NT_coord,scatter_coord,k,Rician_factor);
        [H_Downlink,H_LoS,H_NLoS,NLoS_steering] = GenChannelDL(N_UE,N,L,NR_coord,NT_coord,scatter_coord,k,Rician_factor);
        H_DL_vec=reshape(H_Downlink,[N*N_UE,1]);
        for sampleNum_idx = 1:numel(sampleNumList)
            sampleNum = sampleNumList(sampleNum_idx);


            % sampleNum = sampleNumList(1);
            y = zeros(sampleNum,1);
            Phi = zeros(sampleNum,N*N_UE);
            for m = 1:sampleNum
                F_RF = exp(1j*rand(N,4)*2*pi);
                F_BB = (randn(4,1)+randn(4,1)*1j)/sqrt(2);
                f = F_RF*F_BB;
                f_norm = f/norm(f);
                W_RF = exp(1j*rand(N_UE,1)*2*pi);
                W_BB = 1;
                w = W_RF*W_BB;
                w_norm = w'/norm(w);
                y(m) = w_norm*(H_Downlink*f_norm); % ignore noise in the near-field
                Phi(m,:) = kron(f_norm.',w_norm);
            end
            %         distance = sqrt(UE.y^2+(NT_coord.x-UE.x).^2);
            % y=Phi*H_DL_vec;
            %         norm(y-Phi*H_DL_vec)

            % DFT
            [res_dft,~] = Successive_OMP(y,Phi,Psi_DFT,epsilon,eye(sampleNum,sampleNum),OMP_iternum);
            rec_h_omp_dft = Psi_DFT*res_dft;
            % for iter_ = 1:OMP_iternum
            res_iter = norm(rec_h_omp_dft(:,end)-H_DL_vec)^2/norm(H_DL_vec)^2;
            if isnan(res_iter)
                NMSE_DFT(beta_i,imp,sampleNum_idx) = 0;
            else
                NMSE_DFT(beta_i,imp,sampleNum_idx) = res_iter;
            end
            % end
            % POL
            [ res_pol,~ ] = Successive_OMP( y,Phi,Psi_POL,epsilon,eye(sampleNum,sampleNum),OMP_iternum);
            rec_h_omp_polar = Psi_POL*res_pol;
            % for iter_ = 1:OMP_iternum
            res_iter = norm(rec_h_omp_polar(:,end)-H_DL_vec)^2/norm(H_DL_vec)^2;
            if isnan(res_iter)
                NMSE_POLAR(beta_i,imp,sampleNum_idx) = 0;
            else
                NMSE_POLAR(beta_i,imp,sampleNum_idx) = res_iter;
            end
            % end

            % PROLATE

            err = error_HOLO(imp);
            % we lost error distribution here.

            % Since we have very small x-grid, the error on x-axis is most
            % likely to be small here, div 100 is an approximation to true
            % error distribution.

            x_err = (rand())*err/100;
            y_err = sqrt(err^2-x_err^2);
            x_hat = UE.x+(-1)^randi(2)*x_err;
            y_hat = UE.y+(-1)^randi(2)*y_err;

            distance_hat = sqrt(y_hat^2+(NT_coord.x-x_hat).^2);
            Psi_DPSS_c = Psi_DPSS.'.*exp(-1j*k.*distance_hat);
            Psi_DPSS_c_multiple = kron(conj(Psi_DPSS_c),aR_DPSS);
            for i = 1:L-1
                NLoS_steering(:,i) = NLoS_steering(:,i)/norm(NLoS_steering(:,i),'fro');
            end
            [ res_sph,~ ] = Successive_OMP( y,Phi,[NLoS_steering,Psi_DPSS_c_multiple'],epsilon,eye(sampleNum,sampleNum),OMP_iternum);
            rec_h_omp_sph = [NLoS_steering,Psi_DPSS_c_multiple']*res_sph;
            % for iter_ = 1:OMP_iternum
            res_iter = norm(rec_h_omp_sph(:,end)-H_DL_vec)^2/norm(H_DL_vec)^2;
            if isnan(res_iter)
                NMSE_DPSS(beta_i,imp,sampleNum_idx) = 0.001;
            else
                NMSE_DPSS(beta_i,imp,sampleNum_idx) = res_iter;
            end
            % end
        end
    end
    toc
end
%% RES Analysis
% index_ = 3; % mu = 0.25
% NMSE_DPSS(:,:,1) = NMSE_POLAR(:,:,1);
NMSE_DFT_dB = 10*log10(mean(squeeze(NMSE_DFT),2));
NMSE_POL_dB = 10*log10(mean(squeeze(NMSE_POLAR),2));
NMSE_DPSS_dB = 10*log10(mean(squeeze(NMSE_DPSS),2));
% index = 3;
% INDEX = 1; X RANGE 1~15, Y RANGE -15~0;
% INDEX = 2; X RANGE 1~25, Y RANGE -25~0;
% INDEX = 3; X RANGE 1~25, Y RANGE -40~0;
% err_DFT = std(squeeze(10*log10(NMSE_DFT(index,:,:))),1);
% err_POL = std(squeeze(10*log10(NMSE_POLAR(index,:,:))),1);
% err_SPH = std(squeeze(10*log10(NMSE_SPH(index,:,:))),1);
%%
x = 1:OMP_iternum;
x_sub = 1:2:OMP_iternum;
% errorbar(x,squeeze(NMSE_DFT_dB(index,:,:)),err_DFT,'-s',"LineWidth",1.25,'MarkerSize',8)
% hold on
% errorbar(x,squeeze(NMSE_POL_dB(index,:,:)),err_POL,'-s',"LineWidth",1.25,'MarkerSize',8)
% errorbar(x,squeeze(NMSE_SPH_dB(index,:,:)),err_SPH,'-s',"LineWidth",1.25,'MarkerSize',8)
figure
plot(sampleRatio,squeeze(NMSE_DFT_dB(1,:,:)),'-o',"LineWidth",1.5,'MarkerSize',8)
hold on
plot(sampleRatio,squeeze(NMSE_DFT_dB(2,:,:)),'--o',"LineWidth",1.5,'MarkerSize',8)
plot(sampleRatio,squeeze(NMSE_DFT_dB(3,:,:)),':o',"LineWidth",1.5,'MarkerSize',8)


plot(sampleRatio,squeeze(NMSE_POL_dB(1,:,:)),'-v',"LineWidth",1.5,'MarkerSize',8)
plot(sampleRatio,squeeze(NMSE_POL_dB(2,:,:)),'--v',"LineWidth",1.5,'MarkerSize',8)
plot(sampleRatio,squeeze(NMSE_POL_dB(3,:,:)),':v',"LineWidth",1.5,'MarkerSize',8)

plot(sampleRatio,squeeze(NMSE_DPSS_dB(1,:,:)),'-s',"LineWidth",1.5,'MarkerSize',8)
plot(sampleRatio,squeeze(NMSE_DPSS_dB(2,:,:)),'--s',"LineWidth",1.5,'MarkerSize',8)
plot(sampleRatio,squeeze(NMSE_DPSS_dB(3,:,:)),':s',"LineWidth",1.5,'MarkerSize',8)

grid on
legend('DFT Codebook, \beta=1','DFT Codebook, \beta=2','DFT Codebook, \beta=3','Spherical Codebook, \beta=1','Spherical Codebook, \beta=2','Spherical Codebook, \beta=3','DPSS Codebook, \beta=1','DPSS Codebook, \beta=2','DPSS Codebook, \beta=3')
xlabel('Number of OMP Iterations (I)')
ylabel('NMSE[dB]')
axis([0.01 0.6 -25 0])
if ismac
    set(gca,'fontsize',14);
end

rmpath('../utils/')
if implementations>=1000
    save('../results/NMSE_DFT_dB_red_overhead.mat','NMSE_DFT_dB')
    save('../results/NMSE_POL_dB_red_overhead.mat','NMSE_POL_dB')
    save('../results/NMSE_DPSS_dB_red_overhead.mat','NMSE_DPSS_dB')
end

