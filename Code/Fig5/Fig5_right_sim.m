clear all
close all
rng(42)
%%
addpath('../utils/')
aperture_list = linspace(0.5, 5, 10);
N_UE = 4;
fc = 28e9;
fc_list = [28e9, 39e9, 47e9, 60e9];
c = physconst('Lightspeed');

G_distance = 20;
epsilon = 0;
implementations = 2000;
beta = 1;
error_HOLO = zeros(numel(fc_list), numel(aperture_list));
tic

for aperture_idx = 1:numel(aperture_list)
    aperture = aperture_list(aperture_idx);

    for fc_idx = 1:numel(fc_list)
        fc = fc_list(fc_idx)
        lambda = c / fc;
        d = lambda / 2;
        % aperture = (N_BS-1)*d;
        N_BS = round(aperture / d + 1);

        if mod(N_BS, 2) ~= 0
            N_BS = N_BS + 1;
        end

        k = 2 * pi / lambda;
        array_x = linspace(-aperture / 2, aperture / 2, N_BS);
        array_y = zeros(1, N_BS);
        BS.x = array_x;
        BS.y = array_y;
        UE_range_x = 5;
        UE_range_y_min = 6;
        UE_range_y_max = 20;
        apperture_UE = (N_UE - 1) * d;

        for imple = 1:implementations
            x_s = 2 * (rand() - 0.5) * UE_range_x;
            y_s = UE_range_y_min + (UE_range_y_max - UE_range_y_min) * rand();
            UE.x = x_s + linspace(-apperture_UE / 2, apperture_UE / 2, N_UE);
            UE.y = y_s + zeros(1, N_UE);
            L_NLoS = 3;
            L = L_NLoS + 1;
            Rician_factor = 1;
            recv_pattern = 0;

            for i = 1:numel(UE.x)
                recv_pattern = recv_pattern + cal_recv_pattern(array_x, array_y, UE.x(i), UE.y(i), k, 1) / sqrt(N_UE);
            end

            ref_pattern = 10 * exp(1j * rand(N_BS, 1) * 2 * pi);
            hologram = abs(recv_pattern.' + ref_pattern) .^ 2 - abs(ref_pattern) .^ 2 - abs(recv_pattern.') .^ 2;
            y_grid = linspace(1, 20, G_distance * beta);
            %% fast impl, fft

            x_size = 0;

            while true
                x_size = x_size + 1;

                if x_size * aperture >= UE_range_x
                    %                     x_size
                    x_size = x_size + 1;
                    % if mod(x_size,2)~=0
                    %     x_size = x_size+1;
                    % end
                    break
                end

            end

            HR = ref_pattern .* hologram;
            HR_padded = [zeros((x_size - 0.5) * N_BS, 1); HR; zeros((x_size - 0.5) * N_BS, 1)];
            % HR_padded = [zeros((x_size-0.5)*N_BS,1);HR;zeros((x_size-0.5)*N_BS,1)];
            x_recon_plane = linspace(-aperture * x_size, aperture * x_size, 2 * x_size * N_BS); %
            rec_grid_fft = zeros(numel(x_recon_plane), numel(y_grid));

            for y_idx = 1:numel(y_grid)
                y_s_trg = y_grid(y_idx);
                projection = exp(1j * k * sqrt((x_recon_plane) .^ 2 + y_s_trg ^ 2)); %./sqrt((x_recon_plane).^2+y_s_trg^2)*y_s_trg;
                rad_pattern = fftshift(ifft((fft(HR_padded)) .* (fft(projection.')))) / N_BS;
                rec_grid_fft(:, y_idx) = rad_pattern;
            end

            [~, max_idx] = max(abs(rec_grid_fft(:)));
            [max_row, max_col] = ind2sub(size(rec_grid_fft), max_idx);
            xhat = x_recon_plane(max_row); %+1.36512637125/2;
            yhat = y_grid(max_col);
            error = sqrt((xhat - x_s) ^ 2 + (yhat - y_s) ^ 2);
            error_HOLO(fc_idx, aperture_idx) = error_HOLO(fc_idx, aperture_idx) + error ^ 2;
        end

    end

    toc
end
rmpath('../utils/')
toc
error_HOLO_avg = sqrt(error_HOLO / implementations)
save(['error_HOLO_' num2str(implementations) '_' num2str(beta) '.mat'], 'error_HOLO')
