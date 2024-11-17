clear all
close all
clc

%% the effect of blockage (non-stationary)
%% the effect of wideband bp

addpath("../../utils/")
fc = 28e9;
c = physconst('lightspeed');
lambda = c/fc;
d = lambda/2;
N = 256;
N_UE = 4;
aperture = (N-1)*d;
k = 2*pi/lambda;

BS.x = linspace(-aperture/2,aperture/2,N);
BS.y = zeros(1,N);

UE_range_x = 5;
UE_range_y_min = 6;
UE_range_y_max = 20;
apperture_UE = (N_UE-1)*d;

G_distance = 20;
beta = 8;

x_s = 2*(rand()-0.5)*UE_range_x;
y_s = UE_range_y_min+(UE_range_y_max-UE_range_y_min)*rand();
y_s = y_s/2;
UE.x = x_s + linspace(-apperture_UE/2,apperture_UE/2,N_UE);
UE.y = y_s + zeros(1,N_UE);
L_NLoS = 3;
L = L_NLoS + 1;
Rician_factor=1;
carrier_spacing = 240e3;
N_subc = 100;
k_min = k;
k_max = 2*pi*(fc+(N_subc-1)*carrier_spacing)/c;
k_list = linspace(k_min,k_max,N_subc);
recv_pattern = zeros(numel(k_list),N);

for k_idx = 1:numel(k_list)
    kk = k_list(k_idx);
    for i =1:numel(UE.x)
        recv_pattern(k_idx,:) = recv_pattern(k_idx,:)+cal_recv_pattern(BS.x,BS.y,UE.x(i),UE.y(i),kk,1)/sqrt(N_UE);
    end
end
x_size=0;
while true
    x_size = x_size+1;
    if x_size * aperture >= UE_range_x
        x_size;
        break
    end
end

ref_pattern = 10*exp(1j*rand(N,1)*2*pi);
hologram = abs(recv_pattern(1,:).'+ref_pattern).^2-abs(ref_pattern).^2-abs(recv_pattern(1,:).').^2;
y_grid = linspace(1,20,G_distance*beta);


HR = ref_pattern.*hologram;
HR_padded = [zeros((x_size-0.5)*N,1);HR;zeros((x_size-0.5)*N,1)];
% HR_padded = [zeros((x_size-0.5)*N,1);HR;zeros((x_size-0.5)*N,1)];
x_recon_plane = linspace(-aperture*x_size,aperture*x_size,2*x_size*N);%
rec_grid_fft = zeros(numel(x_recon_plane),numel(y_grid));
for y_idx = 1:numel(y_grid)
    y_s_trg = y_grid(y_idx);
    projection = exp(1j*k*sqrt((x_recon_plane).^2+y_s_trg^2));%./sqrt((x_recon_plane).^2+y_s_trg^2)*y_s_trg;
    rad_pattern = fftshift(ifft((fft(HR_padded)).*(fft(projection.'))))/N;
    rec_grid_fft(:,y_idx) = rad_pattern;
end
block_frac_list = [1/8,1/4,1/2];
block_start = [0,0,0];
rec_grid_fft_block = zeros(numel(x_recon_plane),numel(y_grid),numel(block_frac_list));
for ii = 1:numel(block_frac_list)
    block_frac = block_frac_list(ii);
    N_block = floor(block_frac*N);
    block_start(ii) = round(randi(N-N_block-1));
    hologram_m = hologram;
    hologram_m(block_start(ii):block_start(ii)+N_block) = 0;
    HR = ref_pattern.*hologram_m;
    HR_padded = [zeros((x_size-0.5)*N,1);HR;zeros((x_size-0.5)*N,1)];
    for y_idx = 1:numel(y_grid)
        y_s_trg = y_grid(y_idx);
        projection = exp(1j*k*sqrt((x_recon_plane).^2+y_s_trg^2));%./sqrt((x_recon_plane).^2+y_s_trg^2)*y_s_trg;
        rad_pattern = fftshift(ifft((fft(HR_padded)).*(fft(projection.'))))/N;
        rec_grid_fft_block(:,y_idx,ii) = rad_pattern;
    end
end

figure(1)
surf(x_recon_plane,y_grid,abs(rec_grid_fft)')
shading interp
xlabel('$x$ (m)','Interpreter','Latex');
ylabel('$y$ (m)','Interpreter','Latex');
title('No Blockage','Interpreter','Latex')
set(gca,'fontsize',18);
set(gca,'fontname','Times New Roman')
view([0,0,1])
axis('tight')
box on

figure
surf(x_recon_plane,y_grid,squeeze(abs(rec_grid_fft_block(:,:,1)))')
shading interp
xlabel('$x$ (m)','Interpreter','Latex');
ylabel('$y$ (m)','Interpreter','Latex');
title('$12.5$\% Blockage','Interpreter','Latex')
set(gca,'fontsize',18);
set(gca,'fontname','Times New Roman')
view([0,0,1])
axis('tight')
box on


figure
surf(x_recon_plane,y_grid,squeeze(abs(rec_grid_fft_block(:,:,2)))')
shading interp
xlabel('$x$ (m)','Interpreter','Latex');
ylabel('$y$ (m)','Interpreter','Latex');
title('$25$\% Blockage','Interpreter','Latex')
set(gca,'fontsize',18);
set(gca,'fontname','Times New Roman')
view([0,0,1])
axis('tight')
box on

figure
surf(x_recon_plane,y_grid,squeeze(abs(rec_grid_fft_block(:,:,3)))')
shading interp
xlabel('$x$ (m)','Interpreter','Latex');
ylabel('$y$ (m)','Interpreter','Latex');
title('$50$\% Blockage','Interpreter','Latex')
set(gca,'fontsize',18);
set(gca,'fontname','Times New Roman')
view([0,0,1])
axis('tight')
box on
%% multiple carrier
rec_grid_fft_mc = zeros(numel(x_recon_plane),numel(y_grid),numel(k_list));
for k_idx = 1:numel(k_list)
    hologram_k = abs(recv_pattern(k_idx,:).'+ref_pattern).^2-abs(ref_pattern).^2-abs(recv_pattern(k_idx,:).').^2;
    HR_k = ref_pattern.*hologram_k;
    HR_padded_k = [zeros((x_size-0.5)*N,1);HR_k;zeros((x_size-0.5)*N,1)];
    for y_idx = 1:numel(y_grid)
        y_s_trg = y_grid(y_idx);
        projection = exp(1j*k_list(k_idx)*sqrt((x_recon_plane).^2+y_s_trg^2));%./sqrt((x_recon_plane).^2+y_s_trg^2)*y_s_trg;
        rad_pattern = fftshift(ifft((fft(HR_padded_k)).*(fft(projection.'))))/N;
        rec_grid_fft_mc(:,y_idx,k_idx) = rad_pattern;
    end
end
rec_grid_fft_agg = squeeze(mean(rec_grid_fft_mc,3));

figure
surf(x_recon_plane,y_grid,abs(rec_grid_fft_agg)')
shading interp
xlabel('$x$ (m)','Interpreter','Latex');
ylabel('$y$ (m)','Interpreter','Latex');
title('No Blockage ','Interpreter','Latex')
set(gca,'fontsize',18);
set(gca,'fontname','Times New Roman')
view([0,0,1])
axis('tight')
box on

% block
rec_grid_fft_block = zeros(numel(x_recon_plane),numel(y_grid),numel(block_frac_list), numel(k_list));
for k_idx = 1:numel(k_list)
    hologram_k = abs(recv_pattern(k_idx,:).'+ref_pattern).^2-abs(ref_pattern).^2-abs(recv_pattern(k_idx,:).').^2;
%     HR_k = ref_pattern.*hologram_k;
%     HR_padded_k = [zeros((x_size-0.5)*N,1);HR_k;zeros((x_size-0.5)*N,1)];
    for ii = 1:numel(block_frac_list)
        block_frac = block_frac_list(ii);
        N_block = floor(block_frac*N);
        %     block_start(ii) = round(randi(N-N_block-1));
        hologram_m = hologram_k;
        hologram_m(block_start(ii):block_start(ii)+N_block) = 0;
        HR = ref_pattern.*hologram_m;
        HR_padded = [zeros((x_size-0.5)*N,1);HR;zeros((x_size-0.5)*N,1)];
        for y_idx = 1:numel(y_grid)
            y_s_trg = y_grid(y_idx);
            projection = exp(1j*k_list(k_idx)*sqrt((x_recon_plane).^2+y_s_trg^2));%./sqrt((x_recon_plane).^2+y_s_trg^2)*y_s_trg;
            rad_pattern = fftshift(ifft((fft(HR_padded)).*(fft(projection.'))))/N;
            rec_grid_fft_block(:,y_idx,ii,k_idx) = rad_pattern;
        end
    end
end

rec_grid_fft_block_agg = squeeze(mean(rec_grid_fft_block,4));
figure
surf(x_recon_plane,y_grid,squeeze(abs(rec_grid_fft_block_agg(:,:,1)))')
shading interp
xlabel('$x$ (m)','Interpreter','Latex');
ylabel('$y$ (m)','Interpreter','Latex');
title('$12.5$\% Blockage','Interpreter','Latex')
set(gca,'fontsize',18);
set(gca,'fontname','Times New Roman')
view([0,0,1])
axis('tight')
box on


figure
surf(x_recon_plane,y_grid,squeeze(abs(rec_grid_fft_block_agg(:,:,2)))')
shading interp
xlabel('$x$ (m)','Interpreter','Latex');
ylabel('$y$ (m)','Interpreter','Latex');
title('$25$\% Blockage','Interpreter','Latex')
set(gca,'fontsize',18);
set(gca,'fontname','Times New Roman')
view([0,0,1])
axis('tight')
box on

figure
surf(x_recon_plane,y_grid,squeeze(abs(rec_grid_fft_block_agg(:,:,3)))')
shading interp
xlabel('$x$ (m)','Interpreter','Latex');
ylabel('$y$ (m)','Interpreter','Latex');
title('$50$\% Blockage','Interpreter','Latex')
set(gca,'fontsize',18);
set(gca,'fontname','Times New Roman')
view([0,0,1])
axis('tight')
box on
rmpath("../../utils/")