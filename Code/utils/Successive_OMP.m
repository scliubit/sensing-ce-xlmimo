function [res, iter_num,Support_index_set] = Successive_OMP(y_k_com, Phi, Psi, epsilon, D_w,itermax)

[M_Ns, K] = size(y_k_com);    % size of y_k_com, M_Ns = M*Ns, K is the number of subcarrier
h_v_size = size(Psi,2);
assert(K==1);
res = zeros(h_v_size,itermax);

% Compute the whitened equivalent observation matrix
Upsilon_w = (D_w)'\(Phi*Psi);

% Initialize the residual vectors to the input signal vectors and support estimate, where ..._com contain K columns
y_k_w_com = (D_w)'\y_k_com;
Support_index_set = [];
r_k_com = y_k_w_com;
MSE = 2*epsilon+1;    % Pre-define MSE
iter_num = 0;       % Initialize the number of iteration
% h_v_hat = zeros(h_v_size,K);

while (MSE > epsilon)
    % Distributed Correlation
    c_k_com = Upsilon_w'*r_k_com;

    % Find the maximum projection along the different spaces
    [~, index_p] = max(sum(abs(c_k_com),2));

    % Update the current guess of the common support
    Support_index_set = [Support_index_set; index_p];

    % Project the input signal onto the subspace given by the support using WLS
    xi_hat_com = Upsilon_w(:,Support_index_set)\y_k_w_com;

    % Update residual
    r_k_com = y_k_w_com - Upsilon_w(:,Support_index_set)*xi_hat_com;

    % Compute the current MSE
    MSE = 1/(M_Ns*K)*trace(r_k_com'*r_k_com);

    % Compte the number of iteration
    iter_num = iter_num + 1;
    res(Support_index_set,iter_num) = xi_hat_com;
    if iter_num == itermax
        break
    end
    % if itermax==0
    %     ;
    % else
    %
    % end
end

% assign estimated complex channel gains to the sparse vector
% h_v_hat(Support_index_set,:) = xi_hat_com;

end