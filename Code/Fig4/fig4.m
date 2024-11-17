close all
clear all
clc
rng(18)
addpath("../utils/")
%% demo on how it works

N_BS = 512;
N_UE = 4;
fc = 28e9;
c = physconst('Lightspeed');
lambda = c/fc;
d = lambda/2;
aperture = (N_BS-1)*d;
k = 2*pi/lambda;
array_x = linspace(-aperture/2, aperture/2, N_BS)';
array_y = zeros(N_BS, 1);
UE_range_x = 5;
UE_range_y_min = 6;
UE_range_y_max = 18;
x_s=2*(rand()-0.5)*UE_range_x;
y_s=UE_range_y_min+(UE_range_y_max-UE_range_y_min)*rand();

L_NLoS = 3;
% initialize uplink channel and calculate recv pattern
recv_pattern = cal_recv_pattern(array_x,array_y,x_s,y_s,k,1);
coord_x = [x_s];
coord_y = [y_s];
cnt = 0;
while(1)
    %     UE_range_y_min = 4;
    x_ell = 2*(rand()-0.5)*UE_range_x;
    y_ell = UE_range_y_min+(UE_range_y_max-UE_range_y_min)*rand();
    distances = sqrt((coord_x-x_ell).^2+(coord_y-y_ell).^2);
    if min(distances)<2
        continue
    end
    recv_pattern = recv_pattern+sqrt(1)*cal_recv_pattern(array_x,array_y,x_ell,y_ell,k,1);
    coord_x = [coord_x,x_ell];
    coord_y = [coord_y,y_ell];
    cnt=cnt+1;
    if cnt == L_NLoS
        break
    end
end

G_distance = 20;
Rician_factor = 1;
L = L_NLoS+1;
NT_coord.x = linspace(-aperture/2,aperture/2,N_BS);
NT_coord.y = zeros(1,N_BS);
UE.x=x_s;
UE.y=y_s;
apperture_UE = (N_UE-1)*d;
NR_coord.x = UE.x+linspace(-apperture_UE/2,apperture_UE/2,N_UE);
NR_coord.y = UE.y+zeros(1,N_UE);
scatter_coord.x = coord_x(2:end);
scatter_coord.y = coord_y(2:end);
[H_Downlink,H_LoS,H_NLoS,NLoS_steering] = GenChannelDL(N_UE,N_BS,L,NR_coord,NT_coord,scatter_coord,k,Rician_factor);

ref_pattern = 10*exp(1j*rand(N_BS,1)*2*pi);
hologram = abs(recv_pattern+ref_pattern).^2-abs(ref_pattern).^2-abs(recv_pattern).^2;
y_grid = 1:0.05:23;

%% fast impl, fft
tic
x_size=3;
HR = ref_pattern.*hologram;
HR_padded = [zeros((x_size-0.5)*N_BS,1);HR;zeros((x_size-0.5)*N_BS,1)];
HR_padded = [zeros((x_size-0.5)*N_BS,1);HR;zeros((x_size-0.5)*N_BS,1)];
x_recon_plane = linspace(-aperture*x_size,aperture*x_size,2*x_size*N_BS);%
rec_grid_fft = zeros(numel(x_recon_plane),numel(y_grid));
for y_idx = 1:numel(y_grid)
    y_s_trg = y_grid(y_idx);
    projection = exp(1j*k*sqrt((x_recon_plane).^2+y_s_trg^2));%./sqrt((x_recon_plane).^2+y_s_trg^2)*y_s_trg;
    rad_pattern = fftshift(ifft((fft(HR_padded)).*(fft(projection.'))))/N_BS;
    rec_grid_fft(:,y_idx) = rad_pattern;
end

% figure(1)
figure1 = figure;
% subplot(121)
surf(x_recon_plane,y_grid,abs(rec_grid_fft)')
% colormap(flipud(gray))
hold on
[max_value, max_idx] = max(abs(rec_grid_fft(:)));
[max_row, max_col] = ind2sub(size(rec_grid_fft), max_idx);
x_hat = x_recon_plane(max_row);
y_hat = y_grid(max_col);
% plot3(x_hat,y_hat,max_value,'rv',MarkerSize=15,LineWidth=3)
% title('Time Inversion Near-field Reconstruction ($L=3$)',Interpreter='latex')
shading interp
view([0,0,1])
% view([1,-1,1])

% axis('square')
axis('tight')
axis('equal')
xlabel('$x$ [m]',Interpreter='latex')
ylabel('$y$ [m]',Interpreter='latex')
zlabel('$\vert \tilde{E} ({\bf r})\vert$',Interpreter='latex')
% legend('', 'LoS','',Interpreter='latex')
if ismac
    set(gca,'fontsize',14);
    set(gca,'fontname','Times New Roman')
end
box on
% set(gcf,'renderer','Painters')
% max_value/N_BS*sqrt(x_hat^2+y_hat^2)/100
for i =1:N_BS
    plot3(array_x(i),0,0,'r.')
end
toc

%% detection
thres = max_value/2;
pks_list = [];
locs_list = [];
y_list = [];
w_list = [];
for i=1:numel(y_grid)
    if max(abs(rec_grid_fft(:,i)))<thres
        continue
    else
        [pks,locs,w,p] = findpeaks(abs(rec_grid_fft(:,i)),'MinPeakHeight',thres);
        locs_list = [locs_list;locs];
        y_list = [y_list;zeros(numel(locs),1)+i];
        pks_list = [pks_list;pks];
        w_list = [w_list;w];
    end
end
x_coord = x_recon_plane(locs_list);
y_coord = y_grid(y_list);
X = [x_coord',y_coord'];
[idx,C] = kmeans(X,L_NLoS+1);
for i=1:4
    idx_i = find(idx==i);
    pks_i = pks_list(idx_i);
    w_i = w_list(idx_i);
    [~,idx_max] = max(pks_i);
    [~,idx_min] = min(w_i);
    idx_i_max = idx_i(idx_max);
    idx_i_min = idx_i(idx_min);
    plot3(x_coord(idx_i_max),y_coord(idx_i_max),pks_i(idx_max),'ro',MarkerSize=12,LineWidth=3)
    plot3([x_coord(idx_i_max),x_coord(idx_i_max)],[y_coord(idx_i_max),y_coord(idx_i_max)],[0,pks_i(idx_max)],'r--',LineWidth=2)
end

figure2 = figure;
copyobj(get(figure1, 'Children'), figure2);

axis('tight')
% axis('equal')
xlabel('$x$ [m]',Interpreter='latex')
ylabel('$y$ [m]',Interpreter='latex')
zlabel('$\vert \tilde{E} ({\bf r})\vert$',Interpreter='latex')
% legend('', 'LoS','',Interpreter='latex')
if ismac
    set(gca,'fontsize',14);
    set(gca,'fontname','Times New Roman')
end
box on
view([1,1,1])
% title('Time Inversion Near-field Reconstruction ($L=3$)',Interpreter='latex')
rmpath("../utils/")
