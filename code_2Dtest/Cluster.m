clear;
load('fail_smp.mat');

smp_num = size(fail_smp,1);
cluster_num = 2;
D = size(fail_smp, 2) - 1;

w_smp = sampleWeight(fail_smp, sqrt(D/2));

R = fail_smp(:,end);

fail_smp_o = fail_smp;
fail_smp = bsxfun(@rdivide,fail_smp(:,1:end-1),R);

fprintf('Clustering...\n');
[idx,C]=Wkmeans(fail_smp, cluster_num, w_smp);

Cluster_cnt = zeros(cluster_num,1);
for i = 1:cluster_num
    Cluster_cnt(i) = length(find(idx==i));
end
fprintf('Find %d clusters. Cluster member count: max = %d, min = %d.\n',...
    cluster_num,max(Cluster_cnt), min(Cluster_cnt));


min_norm_points = zeros(cluster_num, D);
fprintf('Locating min norm points...\n');
for i = 1:cluster_num
    fprintf('**********************************\n');
    fprintf('calculating %d in %d min-norm point...\n', i, cluster_num);
    C_k = C(i,:);
    samples_k = fail_smp(idx==i,:);
    R_k = min(R(idx==i));
    d_k = max(Cluster_norm(C_k, samples_k));
    [R_res, sim_num] = min_norm(R_k, C_k, d_k);
    sim_times = sim_times + sim_num;
    C(i,:) = R_res*C(i,:);
    fprintf('R_res = %f.\n',R_res);
end

fprintf('Accumulated simulation times = %d \n', sim_times);
disp('Exit cluster phase.');

save('cluster_res.mat','C','idx');

plot(fail_smp_o(:,1),fail_smp_o(:,2),'*','MarkerSize', 1, 'Color', 'r');
hold;
plot(C(:,1),C(:,2),'o','MarkerSize',10, 'MarkerFaceColor','y');
xlim([-6 6]);
ylim([-6 6]);