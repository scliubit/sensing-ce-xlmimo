clear all
close all
clc
%%

fc = 28e9;
c = physconst('Lightspeed');
lambda = c/fc;
N = 256;
N_UE = 4;
d = lambda/2;
aperture = (N-1)*d; % ~1m
k = 2*pi/lambda;
bandwidth = 9e6;
L=5;
Rician_factor = 1;%~10dB
NT_coord.x = linspace(-aperture/2,aperture/2,N);
NT_coord.y = zeros(1,N);

% UE.x=2*(rand()-0.5)*UE_range_x;
% UE.y=UE_range_y_min+(UE_range_y_max-UE_range_y_min)*rand();
apperture_UE = (N_UE-1)*d;
x_grid = 100;
y_grid = 100;
recv_region_x = linspace(-5,5,x_grid);
recv_region_y = linspace(2,25,y_grid);
dS = (10/x_grid)*(23/y_grid);
recv_region = zeros(x_grid,y_grid);
PT = 0.2; % 23 dBm
if false
    for x_idx = 1:x_grid
        xx = recv_region_x(x_idx);
        for y_idx = 1:y_grid
            yy = recv_region_y(y_idx);
            for ii = 1:N
                d_ii = sqrt((xx-NT_coord.x(ii))^2+(yy-NT_coord.y(ii))^2);
                recv_region(x_idx,y_idx) = recv_region(x_idx,y_idx) + exp(-1j*k*d_ii)/d_ii*PT/N;
            end
        end
    end
end
recv_pwr = abs(recv_region).^2*dS;

%% uplink
random_phase = exp(1j*2*pi*rand(N,1));
recv_region_ul = zeros(x_grid,y_grid);
for x_idx = 1:x_grid
    xx = recv_region_x(x_idx);
    for y_idx = 1:y_grid
        yy = recv_region_y(y_idx);
        for ii = 1:N
            d_ii = sqrt((xx-NT_coord.x(ii))^2+(yy-NT_coord.y(ii))^2);
            recv_region_ul(x_idx,y_idx) = recv_region_ul(x_idx,y_idx) + random_phase(ii)*exp(-1j*k*d_ii)/d_ii*PT/N_UE;
        end
    end
end
recv_pwr_ul = abs(recv_region_ul).^2*lambda/2;
figure
surf(recv_region_x,recv_region_y,10*log10(recv_pwr_ul'*1000))
shading flat
colorbar()
view([1,1,1])
xlabel('$x$ (m)','Interpreter','Latex')
ylabel('$y$ (m)','Interpreter','Latex')
title('Received Power from Near-Field Region (dBm)','Interpreter','Latex')
axis("tight")
% axis("square")
% if ismac
%     set(gca,'fontsize',14);
% end
set(gca,'fontsize',14);
set(gca,'fontname','Times New Roman')
box on
grid on
% 10*log10(mean(recv_pwr_ul'*1000,"all"))
figure
plot(sort(reshape(10*log10(recv_pwr_ul'*1000),[x_grid*y_grid,1])),linspace(0,1,x_grid*y_grid),'-',LineWidth=2,Color="#0072BD")
hold on
plot([-70,20],[0.1,0.1],'r--',LineWidth=2.5)
% histogram()
xlabel('Power (dBm)','Interpreter','Latex')
ylabel('CDF','Interpreter','Latex')
% title('Received Power from Near-Field Region (dBm)','Interpreter','Latex')
set(gca,'fontsize',14);
set(gca,'fontname','Times New Roman')
box on
grid on
return