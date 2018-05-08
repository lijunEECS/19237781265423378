%% only use the vth0 of MP2 and MN5 as variation parameter
% Mean and Sigma are in row-wise
%nmos_vth = 0.50308;
%pmos_vth = -0.4606;
%nmos_sigma = nmos_vth*0.05;
%pmos_sigma = abs(pmos_vth*0.05);
p_mean=[2.7e-9 5.1e-9 1.8e-8 -0.39601 8.80736e-3 -0.15];
n_mean=[2.37e-9 5.8e-9 1.7e-8 0.328977 0.026049 -0.154];
p_sigma=[3.376e-20 4.277e-21 5.687e-20 1.150e-2 4.196e-5 1.797e-3];
n_sigma=[3.602e-22 4.681e-20 1.156e-19 1.094e-2 5.942e-6 1.367e-2];
% k1=size(p_mean,2);
% k2=size(n_mean,2);
% mean_vals=zeros(1,3*(k1+k2)*10);
% sigma_vals=zeros(1,3*(k1+k2)*10);
% mean_vals=zeros(60,k1);
% sigma_vals=zeros(60,k1);
mean_vals=zeros(1,360);
sigma_vals=zeros(1,360);
% for i=1:10
%     mean_vals(3*(k1+k2)*(i-1)+1:3*(k1+k2)*i)=[p_mean p_mean p_mean n_mean n_mean n_mean];
%     sigma_vals(3*(k1+k2)*(i-1)+1:3*(k1+k2)*i)=[p_sigma p_sigma p_sigma n_sigma n_sigma n_sigma];
% end
for i=1:10
    mean_vals(:,36*(i-1)+1:36*(i-1)+18)=[p_mean p_mean p_mean];
    mean_vals(:,36*(i-1)+19:36*(i-1)+36)=[n_mean n_mean n_mean];
    sigma_vals(:,36*(i-1)+1:36*(i-1)+18)=[p_sigma p_sigma p_sigma];
    sigma_vals(:,36*(i-1)+19:36*(i-1)+36)=[n_sigma n_sigma n_sigma];
end
