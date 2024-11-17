function [Psi_DFT,Psi_POL,Psi_DPSS,index2coor] = dict_design_red(N_UE,N,aperture,k,z_min,z_max,x_max,beta)
z_max = z_max+2;
x_range = linspace(-aperture/2,aperture/2,N);
aperture_UE = aperture/N*N_UE;
% lambda = 2/pi/k;
% DFT DICT
if N_UE>1
    aR_DFT = dftmtx(N_UE*beta)/sqrt(N_UE*beta);
    aR_DFT = aR_DFT(1:N_UE,:);
    aT_DFT = dftmtx(N*beta)/sqrt(N*beta);
    aT_DFT = aT_DFT(1:N,:);
    Psi_DFT = kron((aT_DFT),aR_DFT);
else
    Psi_DFT = dftmtx(N*beta)/sqrt(N*beta);
    Psi_DFT = Psi_DFT(1:N,:);
end
% WAVE NUMBER

% i_BS_min = ceil(-aperture/lambda);
% i_BS_max = floor(aperture/lambda);
% i_BS = i_BS_min:1/beta:i_BS_max;
% if N_UE>1
%     i_UE_min = ceil(-aperture_UE/lambda);
%     i_UE_max = floor(aperture_UE/lambda);
%     i_UE = i_UE_min:1/beta:i_UE_max;
%     r = linspace(-aperture,aperture,N)';
%     Psi_BS = exp(1j*2*pi*i_BS*r);
%     rue = linspace(-aperture,aperture,N)';
% else
% end

% POL DICT
if true
if N_UE>1
    r_inv = linspace(1/z_max,1/z_min,round(sqrt(N*beta)));%0:0.01:1;
    r = 1./r_inv;
    x = linspace(-x_max,x_max,round(sqrt(N*beta)));
    aT_POL = zeros(numel(x_range),numel(r_inv)*numel(x));
    cnt=0;
    index2coor=zeros(numel(r_inv)*numel(x),2);
    for xx =1:numel(x)
        for rr = 1:numel(r)
            cnt=cnt+1;
            %         DICT_Layer1(:,cnt)=exp(-1j*k.*(r(rr)+(x_range.'-x(xx)).^2/(r(rr))));
            aT_POL(:,cnt)=exp(-1j*k.*sqrt(r(rr)^2+(x(xx)-x_range.').^2));%./sqrt(r(rr)^2+(x(xx)-x_range.').^2);
            index2coor(cnt,1) = r(rr);
            index2coor(cnt,2) = x(xx);
        end
    end
    index2coor = kron(index2coor,[1;1;1;1]);
    % ==============================
    x_range_R = linspace(-aperture_UE/2,aperture_UE/2,N_UE);
    r_inv = linspace(1/z_max,1/z_min,floor(sqrt(N_UE*beta)));%0:0.01:1;
    r = 1./r_inv;
    x = linspace(-x_max,x_max,ceil(sqrt(N_UE*beta)));
    aR_POL = zeros(numel(x_range_R),numel(r_inv)*numel(x));
    cnt=0;
    index2coor2=zeros(numel(r_inv),2);
    for xx =1:numel(x)
        for rr = 1:numel(r)
            cnt=cnt+1;
            % r_true = sqrt(rr^2+xx^2);
            %         DICT_Layer1(:,cnt)=exp(-1j*k.*(r(rr)+(x_range.'-x(xx)).^2/(r(rr))));
            aR_POL(:,cnt)=exp(-1j*k.*sqrt(r(rr)^2+(x(xx)-x_range_R.').^2));%./sqrt(r(rr)^2+(x(xx)-x_range.').^2);
            index2coor2(cnt,1) = r(rr);
            index2coor2(cnt,2) = x(xx);
        end
    end
    Psi_POL = kron(aT_POL,aR_POL);
    Psi_POL = kron(aT_POL,aR_DFT);
else
    r_inv = linspace(1/z_max,1/z_min,floor(sqrt(N*beta)));%0:0.01:1;
    r = 1./r_inv;
    x = linspace(-x_max,x_max,ceil(sqrt(N*beta)));
    Psi_POL = zeros(numel(x_range),numel(r_inv)*numel(x));
    cnt=0;
    index2coor=zeros(numel(r_inv),2);
    for xx =1:numel(x)
        for rr = 1:numel(r)
            cnt=cnt+1;
            %         DICT_Layer1(:,cnt)=exp(-1j*k.*(r(rr)+(x_range.'-x(xx)).^2/(r(rr))));
            Psi_POL(:,cnt)=exp(-1j*k.*sqrt(r(rr)^2+(x(xx)-x_range.').^2));%./sqrt(r(rr)^2+(x(xx)-x_range.').^2);
            index2coor(cnt,1) = r(rr);
            index2coor(cnt,2) = x(xx);
        end
    end
end
end
% POL 2
if false
if N_UE>1
    theta_ = linspace(-0.93,0.93,round(sqrt(N)*beta));
    r_inv = linspace(1/z_max,1/z_min,round(sqrt(N)*beta));%0:0.01:1;
    r = 1./r_inv;
%     x = linspace(-x_max,x_max,round(sqrt(N)*beta));
    aT_POL = zeros(N,numel(r_inv)*numel(theta_));
    cnt=0;
    index2coor=zeros(numel(r_inv)*numel(theta_),2);
    for xx =1:numel(theta_)
        for rr = 1:numel(r)
            cnt=cnt+1;
            [x_,y_] = pol2cart(acos(theta_(xx)),r(rr)*(1-theta_(xx)^2));
            %         DICT_Layer1(:,cnt)=exp(-1j*k.*(r(rr)+(x_range.'-x(xx)).^2/(r(rr))));
            aT_POL(:,cnt)=exp(-1j*k.*sqrt(y_^2+(x_-x_range.').^2));%./sqrt(r(rr)^2+(x(xx)-x_range.').^2);
            index2coor(cnt,1) = y_;
            index2coor(cnt,2) = x_;
        end
    end
    index2coor = kron(index2coor,[1;1;1;1]);
    % ==============================
    theta_R = linspace(-0.93,0.93,ceil(sqrt(N_UE*beta)));
    x_range_R = linspace(-aperture_UE/2,aperture_UE/2,N_UE);
    r_inv = linspace(1/z_max,1/z_min,ceil(sqrt(N_UE*beta)));%0:0.01:1;
    r = 1./r_inv;
%     x = linspace(-x_max,x_max,ceil(sqrt(N_UE)));
    aR_POL = zeros(N_UE,numel(r_inv)*numel(theta_R));
    cnt=0;
    index2coor2=zeros(numel(r_inv),2);
    for xx =1:numel(theta_R)
        for rr = 1:numel(r)
            cnt=cnt+1;
            [x_,y_] = pol2cart(acos(theta_R(xx)),r(rr)*(1-theta_R(xx)^2));
            %         DICT_Layer1(:,cnt)=exp(-1j*k.*(r(rr)+(x_range.'-x(xx)).^2/(r(rr))));
            aR_POL(:,cnt)=exp(-1j*k.*sqrt(y_^2+(x_-x_range_R.').^2));%./sqrt(r(rr)^2+(x(xx)-x_range.').^2);
            index2coor2(cnt,1) = y_;
            index2coor2(cnt,2) = x_;
        end
    end
    Psi_POL = kron(aT_POL,aR_POL);
    Psi_POL = kron(aT_POL,eye(N_UE,N_UE));
else
    theta_ = linspace(-0.93,0.93,ceil(N*beta/2));
    r_inv = linspace(1/z_max,1/z_min,10);%0:0.01:1;
    r = 1./r_inv;
%     x = linspace(-x_max,x_max,ceil(sqrt(N*beta)));
    Psi_POL = zeros(N,numel(r_inv)*numel(theta_));
    cnt=0;
    index2coor=zeros(numel(r_inv),2);
    for xx =1:numel(theta_)
        for rr = 1:numel(r)
            cnt=cnt+1;
            [x_,y_] = pol2cart(acos(theta_(xx)),r(rr)*(1-theta_(xx)^2));
            %         DICT_Layer1(:,cnt)=exp(-1j*k.*(r(rr)+(x_range.'-x(xx)).^2/(r(rr))));
            Psi_POL(:,cnt)=exp(1j*k.*sqrt(y_^2+(x_-x_range.').^2));%./sqrt(r(rr)^2+(x(xx)-x_range.').^2);
            index2coor(cnt,1) = y_;
            index2coor(cnt,2) = x_;
        end
    end
end

end

% Prolate DICT


if N_UE>1
    time_halfbandwidth = 0;%k*aperture/4/pi/z_max;
    time_halfbandwidth_R = 1;%k*aperture_UE/4/pi/z_max;
    % num_seq = 4;
    [aT_DPSS,~] = dpss(N,time_halfbandwidth,N,"trace");
    [aR_DPSS,~] = dpss(N_UE,time_halfbandwidth_R,N_UE,"trace");
    Psi_DPSS = aT_DPSS;%kron(conj(aT_DPSS),aR_DPSS);
else
    time_halfbandwidth = 4;%k*aperture/4/pi/z_max;
    % num_seq = 4;
    [Psi_DPSS,~] = dpss(N,time_halfbandwidth,N,"trace");
end

end

