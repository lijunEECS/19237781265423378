function [MCpfail, MCfom, sample_n, new_smp, sim_times] = mis_main(fail_smp, C, idx, sim_times)

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
stop_fom = 0.1;
MCtotal_error_counter = 0;
MCtotal_weight_sum = 0;
sample_n=0;
iter=1;
sample_unit = 2000;
threshold = 139.5;
new_smp = [];
max_iter = 50;


disp('**********************************************');
disp('Deploy Mixture Importance Sampling Monte Carlo Simulation...');


while(iter<max_iter)
    
    [samples, w_smp] = generateMISSamples(C, alpha, beta, sample_unit);
    
    MCresults = isFailure(samples, threshold);
    MCerror_counter = nnz(MCresults) ;
    MCtotal_error_counter = MCtotal_error_counter + MCerror_counter;
    sample_n = [sample_n, sample_n(end)+sample_unit];
    
    error_idx = find(MCresults);
    MCweight_sum = sum(w_smp(error_idx));
    MCtotal_weight_sum = MCtotal_weight_sum + MCweight_sum;
    Rs = sqrt(sum(samples.^2,2));
    Stmp = [samples(error_idx,:) Rs(error_idx)];
    new_smp = [new_smp; Stmp];
    [N,D] = size(new_smp);
    new_smp = sortrows(new_smp,D);
    N = min(N, 1000);
    new_smp = new_smp(1:N,:);
    
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
            return;
        end
    end
    if(MCfom(end)<0.20)
        iter = 1;
    end
end


end

