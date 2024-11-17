function [recv_pattern,norm_vec] = cal_recv_pattern(array_x,array_y,x_s,y_s,kappa,dir)
%% back proj
array_dir = [array_x(1)-array_x(end),array_y(1)-array_y(end)];
array_dir = array_dir/norm(array_dir);
norm_vec = [-array_dir(2),array_dir(1)];
distance = sqrt((array_x-x_s).^2+(array_y-y_s).^2);
cos_theta = ((array_x-x_s)*norm_vec(1)+(array_y-y_s)*norm_vec(2))./distance;
recv_pattern = exp(-1*dir*1j*kappa*distance)./distance;
end