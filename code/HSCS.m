% seed=20180429;
%
% rng(seed);

D=360;
dr = sqrt(D)/3;
R = dr;
fail_smp = [];
fail_cnt = 0;
sample_unit = 2000;
threshold = 139.5;
sub_threshold = 138.0;
t_threshold=0.10;
p_threshold=0.01;
t=0;
fprintf('Deploy spherical presampling...\n');
sim_times=0;
iter = 0;
find_flag = 0;

while(1)
    r_samples = normrnd(0,1,sample_unit, D);
    r_norm = sqrt(sum(abs(r_samples).^2,2));
    r_normalized_samples_o = bsxfun(@rdivide,r_samples,r_norm);
    
    r_normalized_samples = R*r_normalized_samples_o;
    
    p=2*normcdf(-abs(r_normalized_samples));
    
    t=sum(p(:)<p_threshold)/(sample_unit*D);
    
    sim_times = sim_times+sample_unit;
    fprintf('==============================\n');
    fprintf('Accumulated simulation times = %d \n', sim_times);
    fprintf('t=%f \n', t);
    
    samples = r_normalized_samples;
    
    if(~find_flag)
        MCresults =isFailure(samples, sub_threshold);
        fail_idx = find(MCresults);    
         if(isempty(fail_idx))
            fprintf('R=%f, didnt find failed samples.\n', R);
         else
            dr = sqrt(D)/10;
            find_flag = 1;
            sample_unit = 20000;
            fprintf('R=%f, get into suspicious region.\n', R);
            R = R - dr;
         end
    else
        MCresults =isFailure(samples, threshold);

        fail_idx = find(MCresults);
    
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
    end
    
    if(t>t_threshold && fail_cnt>250)
        disp('Exit presampling...');
        break;
    end
    
    if (t<t_threshold)
        R=R+dr;
    end
    iter = iter +1;
end

save('presample_res.mat', 'fail_smp', 'sim_times','dr');

%% cluster

% clear;
% load('presample_res.mat');

[C,idx, sim_times] = kmeans_main(fail_smp,dr,sim_times);

save('cluster_res.mat','C','idx','sim_times','fail_smp');

%% mix importance sampling

clear;
load('cluster_res.mat');

[MCpfail, MCfom, sample_n, new_smp, sim_times] = mis_main(fail_smp,C,idx, sim_times);
file_name = [num2str(0), 'thMIS_res.mat'];
save(file_name);
    
stop_fom = 0.1;
iter=1;
while(MCfom(end)>stop_fom)
    dr = (new_smp(end,end)-new_smp(1,end))/5;
    [C, idx, sim_times] = kmeans_main(new_smp, dr, sim_times);
    fprintf('Accumulated simulation times = %d \n', sim_times);
    [MCpfail, MCfom, sample_n, new_smp, sim_times] = mis_main(new_smp,C,idx, sim_times);
    file_name = [num2str(iter), 'thMIS_res.mat'];
    save(file_name);
    iter = iter+1;
end

fprintf('Accumulated simulation times = %d.\n', sim_times);
fprintf('Exit mixIS phase.\n');

figure(1)
semilogx(sample_n(2:end),MCpfail, '-*'); title('MCPFail');
xlabel('number of samples');
figure(2)
semilogx(sample_n(2:end),MCfom, '-*'); title('MCfom');
xlabel('number of samples');

save('mixISres.mat');


