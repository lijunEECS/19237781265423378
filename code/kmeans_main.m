function [C, idx, sim_times] = kmeans_main(fail_smp, dr, sim_times)

smp_num = size(fail_smp,1);
D = size(fail_smp, 2) - 1;

w_smp = sampleWeight(fail_smp, dr);

R = fail_smp(:,end);

fail_smp = bsxfun(@rdivide,fail_smp(:,1:end-1),R);

fprintf('Clustering...\n');
max_cluster_num = ceil(sqrt(smp_num));
% max_cluster_num = 60;
idx = zeros(max_cluster_num, smp_num);
C = zeros(max_cluster_num, max_cluster_num, D);
DBIres = zeros(max_cluster_num-1,1);
cluster_num_res = zeros(max_cluster_num-1,1);
for i = 2:max_cluster_num
    tmp_DBIres = zeros(100,1);
    tmp_idx = zeros(100, smp_num);
    tmp_C = zeros(100, max_cluster_num, D);
    tmp_cluster_num = zeros(100,1);
    for j = 1:100
        [tmp_idx(j,:),tmp_C_resize,tmp_cluster_num(j)]=Wkmeans(fail_smp, i, w_smp);
        
        if(size(tmp_C_resize,1)~=tmp_cluster_num(j))
            disp('error');
        end
        
        tmp_C(j,1:tmp_cluster_num(j),:) = tmp_C_resize;
        tmp_DBIres(j) = DBI(fail_smp, squeeze(tmp_idx(j,:)), squeeze(tmp_C(j,1:tmp_cluster_num(j),:)));
    end
    best_idx = find(tmp_DBIres == min(tmp_DBIres));
    best_idx = best_idx(1);
    DBIres(i-1) = tmp_DBIres(best_idx);
    idx(i,:) = tmp_idx(best_idx,:);
    C(i,1:i,:) = tmp_C(best_idx, 1:i, :);
    cluster_num_res(i-1) = tmp_cluster_num(best_idx);
end
    
%     DBIres(i - 1) = DVI(fail_smp, squeeze(idx(i,:)));

cluster_num = cluster_num_res(DBIres==min(DBIres));
cluster_num = cluster_num(end);


idx = idx(cluster_num,:);
C = C(cluster_num, 1:cluster_num, :);
C = squeeze(C);

Cluster_cnt = zeros(cluster_num,1);
for i = 1:cluster_num
    Cluster_cnt(i) = length(find(idx==i));
end
fprintf('Find %d clusters. Cluster member count: max = %d, min = %d.\n',...
    cluster_num,max(Cluster_cnt), min(Cluster_cnt));

fprintf('Locating min norm points...\n');

for i = 1:cluster_num
    fprintf('**********************************\n');
    fprintf('calculating %d in %d min-norm point...\n', i, cluster_num);
    C_k = C(i,:);
    samples_k = fail_smp(idx==i,:);  
    R_k = min(R(idx==i));    
    d_k = max(Cluster_norm(C_k, samples_k));
    if(size(fail_smp,1) <= 5)
        R_res = R_k;
    else
        [R_res, sim_num] = min_norm(R_k, C_k, d_k);
        sim_times = sim_times + sim_num;
    end
    
    %%%test
    test1 = R_res*C(i,:);
    if(isempty(test1))
        disp('error');
    end
    %%%test
    
    C(i,:) = R_res*C(i,:);
    fprintf('R_res = %f.\n',R_res);
end

fprintf('Accumulated simulation times = %d \n', sim_times);
disp('Exit cluster phase.');

end

