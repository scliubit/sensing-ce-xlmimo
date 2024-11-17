clear all
close all
%% MUSIC algorithm takes veeeery long time.

addpath('../utils/')
fc = 28e9;
c = physconst('Lightspeed');
lambda = c/fc;
N = 256;
N_UE = 1;
d = lambda/2;
aperture = (N-1)*d; % ~1m
k = 2*pi/lambda;
bandwidth = 9e6;
L=5;
Rician_factor = 1;
NT_coord.x = linspace(-aperture/2,aperture/2,N);
NT_coord.y = zeros(1,N);
% H_downlink = zeros(1,N);
% UE location range
UE_range_x = 2;
UE_range_y_min = 2;
UE_range_y_max = 20;
implementations = 1000;
G_distance = 20;
G_angle = N*8;
% sampleRatio=[0.25, 0.4, 0.6];
epsilon = 0;
% beta = 1;
for beta = [1,2,3]
    % implementations = 2000, beta = 1;     135.91s
    % implementations = 2000, beta = 2;     404.38s
    % implementations = 2000, beta = 3;     822.24s
    % implementations = 2000, beta = 4;     1289.96s
    error_MUSIC = zeros(1,numel(implementations));
    error_OMP = zeros(1,numel(implementations));
    error_HOLO = zeros(1,numel(implementations));
    [~,Psi_POL,~,index2coor] = dict_design_red(N_UE,N,aperture,k,UE_range_y_min,UE_range_y_max,UE_range_x,beta);
    tic
    for imple = 1:implementations
        if mod(imple,500)==0
            toc
            disp(['[' num2str(imple) ']/[' num2str(implementations) ']'])
        end
        UE.x=2*(rand()-0.5)*UE_range_x;
        UE.y=UE_range_y_min+(UE_range_y_max-UE_range_y_min)*rand();
        apperture_UE = (N_UE-1)*d;
        NR_coord.x = UE.x+linspace(-apperture_UE/2,apperture_UE/2,N_UE);
        NR_coord.y = UE.y+zeros(1,N_UE);
        if L>=2
            theta = linspace(-pi,pi,L-1);
            R = 5;
            scatter_coord.x = R*cos(theta);
            scatter_coord.y = R*sin(theta);
        else
            scatter_coord.x = [];
            scatter_coord.y = [];
        end
        [H_Downlink,H_LoS,H_NLoS,NLoS_steering] = GenChannelDL(N_UE,N,L,NR_coord,NT_coord,scatter_coord,k,Rician_factor);
        %     norm(H_NLoS,'fro')^2/norm(H_LoS,'fro')^2
        % MUSIC
        H_Upink = H_Downlink';

        sampleNum = ceil(0.5*N);
        y_ = zeros(sampleNum,1);
        Phi = zeros(sampleNum,N*N_UE);
        for m = 1:sampleNum
            F_RF = exp(1j*rand(N,1)*2*pi);
            F_BB = (randn(1,1)+randn(1,1)*1j)/sqrt(2);
            f = F_RF*F_BB;
            f_norm = f/norm(f);
            W_RF = exp(1j*rand(N_UE,1)*2*pi);
            W_BB = 1;
            w = W_RF*W_BB;
            w_norm = w'/norm(w);
            y_(m) = w_norm*(H_Downlink*f_norm); % ignore noise in the near-field
            Phi(m,:) = kron(f_norm.',w_norm);
        end
        %% MUSIC
        %         if false
        %         Phi = eye(256,256);
        R = (Phi'*Phi*H_Upink)*(Phi'*Phi*H_Upink)';
        % perform MUSIC
        [~,~,V] = svd(R');
        V = V(:,L+1:end);
        theta = linspace(0.1*pi,0.9*pi,4*N*beta);
        distance_inv = linspace(1/20,1,G_distance*beta);
        distance = 1./distance_inv;
        P_MUSIC = zeros(length(theta),length(distance_inv));
        nT = (0:N-1)';
        for i = 1:length(theta)
            for j = 1:length(distance_inv)
                [x,y] = pol2cart(theta(i),distance(j));
                distance_i = sqrt((NT_coord.x-x).^2+(NT_coord.y-y).^2);
                a = exp(-1j*k*distance_i);
                P_MUSIC(i,j) = abs(1/(a*V*V'*a'));
            end
        end
        [~,index] = max(P_MUSIC(:));
        [angle_index,distance_index] = ind2sub(size(P_MUSIC),index);
        angle_hat = theta(angle_index);
        distance_hat = distance(distance_index);
        [x_hat,y_hat] = pol2cart(angle_hat,distance_hat);
        error = sqrt((x_hat-UE.x)^2+(y_hat-UE.y)^2);
        error_MUSIC(imple) = error;

        %% OMP
        [ h_v_hat,~ ] = OMP_Algorithm_MMV( y_,Phi,conj(Psi_POL),epsilon,eye(sampleNum,sampleNum),1);
        [~,index]= max(abs(h_v_hat));
        y_hat = index2coor(index,1);
        x_hat = index2coor(index,2);
        error = sqrt((x_hat-UE.x)^2+(y_hat-UE.y)^2);
        error_OMP(imple) = error;
        %         end

        % hologram
        recv_pattern = 0;
        for i = 1:N_UE
            recv_pattern = recv_pattern+cal_recv_pattern(NT_coord.x',NT_coord.y',NR_coord.x(i),NR_coord.y(i),k,1)/sqrt(N_UE);
        end
        ref_pattern = 10*exp(1j*rand(N,1)*2*pi);
        ref_pattern = awgn(ref_pattern,10,'measured'); % wont affect results
        hologram = abs(recv_pattern+ref_pattern).^2-abs(ref_pattern).^2-abs(recv_pattern).^2;
        y_grid = linspace(UE_range_y_min,UE_range_y_max,G_distance*beta);
        x_size=4;
        HR = ref_pattern.*hologram;
        x_recon_plane = linspace(-aperture*x_size,aperture*x_size,2*x_size*N);%
        % x_recon_plane = linspace(-aperture*x_size,aperture*x_size,beta*x_size*N);%
        HR_padded = [zeros((x_size-0.5)*N,1);HR;zeros((x_size-0.5)*N,1)];
        rec_grid_fft = zeros(numel(x_recon_plane),numel(y_grid));
        for y_idx = 1:numel(y_grid)
            y_s_trg = y_grid(y_idx);
            projection = exp(1j*k*sqrt((x_recon_plane).^2+y_s_trg^2));%./sqrt((x_recon_plane).^2+y_s_trg^2)*y_s_trg;
            rad_pattern = fftshift(ifft((fft(HR_padded)).*(fft(projection.'))))/N;
            rec_grid_fft(:,y_idx) = rad_pattern;
        end

        [~, max_idx] = max(abs(rec_grid_fft(:)));
        [max_row, max_col] = ind2sub(size(rec_grid_fft), max_idx);
        xhat = x_recon_plane(max_row);%+1.36512637125/2;
        yhat = y_grid(max_col);
        error = sqrt((xhat - UE.x)^2+(yhat - UE.y)^2);
        error_HOLO(imple) = error;
    end
    toc
    %%
    if implementations>=1000
        save(['../results/error_MUSICc_' num2str(implementations) '_' num2str(beta) '.mat'],'error_MUSIC')
        save(['../results/error_OMPc_' num2str(implementations) '_' num2str(beta) '.mat'],'error_OMP')
        save(['../results/error_HOLOc_' num2str(implementations) '_' num2str(beta) '.mat'],'error_HOLO')
    end
end
rmpath('../utils/')
