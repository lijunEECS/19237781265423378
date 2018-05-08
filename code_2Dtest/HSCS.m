function [HSCSres, simulation_times, cluster_num] = HSCS(plotFlag)

% seed=20180429;
%
% rng(seed);

data1 = rand([4000,2]);
data1 = 16*data1;
data1 = data1 - 8;

res = isFailure(data1);
idx1 = find(res);
idx2 = find(~res);

% plot(data1(idx1,1),data1(idx1,2), '*', 'MarkerSize', 5, 'Color','r')
% hold;
% plot(data1(idx2,1),data1(idx2,2), '*', 'MarkerSize', 1, 'Color','b')
% legend('Fail','Accept');
% xlabel('x1');
% ylabel('x2');

D=2;
dr = sqrt(D);
R = dr;
fail_smp = [];
fail_cnt = 0;
sample_unit = 200;
t_threshold=0.05;
p_threshold=2*normcdf(-5);
t=0;
fprintf('Deploy spherical presampling...\n');
sim_times=0;
data=[];
data_fail_idx=[];
iter = 0;
while(1)
    r_samples = normrnd(0,1,sample_unit, D);
    r_norm = sqrt(sum(abs(r_samples).^2,2));
    r_normalized_samples_o = bsxfun(@rdivide,r_samples,r_norm);
    
    r_normalized_samples = R*r_normalized_samples_o;
    
    p=2*normcdf(-abs(r_normalized_samples));
    
    t=sum(p(:)<p_threshold)/(sample_unit*D);
    
    if(t>t_threshold && fail_cnt>500)
        disp('Exit presampling...');
        break;
    end
    
    sim_times = sim_times+sample_unit;
    fprintf('Accumulated simulation times = %d \n', sim_times);
    
    samples = r_normalized_samples;
    data = [data;samples];
    
    MCresults =isFailure(samples);
    
    fail_idx = find(MCresults);
    
    data_fail_idx = [data_fail_idx;MCresults];
    if(isempty(fail_idx))
        fprintf('R=%f, didnt find failed samples.\n', R);
    else
        fail_cnt = fail_cnt + length(fail_idx);
        fprintf('R=%f, fail_cnt=%d.\n', R, fail_cnt);
        w=R*ones(length(fail_idx),1);
        fail_smp = [fail_smp; r_normalized_samples(fail_idx,:) w];
        %         sample_unit = 10*sample_unit_o;
        %         dr=sqrt(sqrt(D/2));
    end
    if (t<t_threshold)
        R=R+dr;
    end
    iter = iter +1;
end

%% cluster

smp_num = size(fail_smp,1);
D = size(fail_smp, 2) - 1;

w_smp = sampleWeight(fail_smp, sqrt(D/2));

R = fail_smp(:,end);

fail_smp = bsxfun(@rdivide,fail_smp(:,1:end-1),R);

fprintf('Clustering...\n');
max_cluster_num = ceil(sqrt(smp_num));
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
    [R_res, sim_num] = min_norm(R_k, C_k, d_k);
    sim_times = sim_times + sim_num;
    
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

if(plotFlag)
    save('cluster_res.mat','C','idx');
    
    figure;
    plot(data1(idx1,1),data1(idx1,2), '*', 'MarkerSize', 5, 'Color','r')
    hold;
    plot(data1(idx2,1),data1(idx2,2), '*', 'MarkerSize', 1, 'Color','b')
    plot(C(:,1),C(:,2),'o','MarkerSize',10, 'MarkerFaceColor','y');
    legend('Fail','Accept','min-norm points');
    xlabel('x1');
    ylabel('x2');
end

%% mix importance sampling

alpha = 10^(-3);
D=2;

w = sampleWeight(fail_smp, sqrt(D/2));
w_sum = sum(w);
cluster_num = size(C,1);
beta = zeros(cluster_num,1);
for i =1:cluster_num
    w_k = w(idx==i);
    beta(i) = sum(w_k)/w_sum;
end

MCpfail = [];
MCfom = [];
stop_fom = 0.05;
MCtotal_error_counter = 0;
MCtotal_weight_sum = 0;
sample_n=0;
iter=0;
sample_unit = 100;
threshold = 139.5;

%C = 0.9*C;

disp('**********************************************');
disp('Deploy Mixture Importance Sampling Monte Carlo Simulation...');


while(1)
    
    [samples, w_smp] = generateMISSamples(C, alpha, beta, sample_unit);
    
    MCresults = isFailure(samples);
    MCerror_counter = nnz(MCresults) ;
    MCtotal_error_counter = MCtotal_error_counter + MCerror_counter;
    sample_n = [sample_n, sample_n(end)+sample_unit];
    
    error_idx = find(MCresults);
    MCweight_sum = sum(w_smp(error_idx));
    MCtotal_weight_sum = MCtotal_weight_sum + MCweight_sum;
    
    iter = iter+1;
    %     if(iter<25)
    %         continue;
    %     end
    
    MCpfail = [MCpfail MCtotal_weight_sum/sample_n(end)];
    MCfom = [MCfom std(MCpfail)/mean(MCpfail)];
    str = sprintf('%d out of %d samples failed(%d/%d), MC failure rate = %e, MC FOM = %e', MCtotal_error_counter, sample_n(end), MCerror_counter, sample_unit, MCpfail(end), MCfom(end));
    disp(str);
    
    sim_times = sim_times + sample_unit;
    
    if(iter > 10)
        if(MCfom(end)<=stop_fom)
            break;
        end
    end
    
end

fprintf('Accumulated simulation times = %d.\n', sim_times);
fprintf('Exit mixIS phase.\n');

if(plotFlag)
    save('fail_smp.mat','fail_smp','sim_times');
    figure;
    plot(data(data_fail_idx>0,1),data(data_fail_idx>0,2),'*','MarkerSize',5,'Color','r');
    hold;
    plot(data(~(data_fail_idx>0),1),data(~(data_fail_idx>0),2),'x','MarkerSize',1,'Color','b');
    legend('Fail','Accept');
    xlabel('x1');
    ylabel('x2');
end

if(plotFlag)
    figure
    semilogx(sample_n(2:end),MCpfail, '-*'); title('MCPFail');
    xlabel('number of samples');
    figure
    semilogx(sample_n(2:end),MCfom, '-*'); title('MCfom');
    xlabel('number of samples');
    
    save('mixISres.mat');
end
HSCSres = MCpfail(end);
simulation_times = sim_times;
end