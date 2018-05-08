function [R_res, sim_num] = min_norm(R, C, d)

    R_max = R;
    R_min = 0;
    R_threshold = R/100;
    batch_size = 1000;
    D = length(C);
    sim_num = 0;
    R_res = R;
    threshold = 139.5;
    if(d<2e-32)
        return;
    end
    R_tmp = R_max;
    while(R_max-R_min>R_threshold) 
        R_tmp = (R_max + R_min)/2;
        samples = [];
        batch_num = 0;
        iter = 0;
        try_size = 10*batch_size;
        while(batch_num < batch_size)
            iter = iter+1;
            sample_try = rand([try_size, D]);
            sample_try = bsxfun(@times, sample_try, 2*C);
            sample_try = bsxfun(@minus, sample_try, C);
            sample_try_norm = sqrt(sum(sample_try.^2,2));
            sample_try = bsxfun(@rdivide, sample_try, sample_try_norm);
            sample_try = R_tmp*sample_try; 

%             sample_try = normrnd(0,1,10*batch_size, D);
%             sample_try = bsxfun(@plus, sample_try, C);
%             sample_try_norm = sqrt(sum(sample_try.^2, 2));
%             sample_try = bsxfun(@rdivide, sample_try, sample_try_norm);
%             sample_try = R_tmp*sample_try;

            dis = Cluster_norm(C, sample_try);
            idx = find(dis<d);
            if(iter>100)
                fprintf('d = %d. Only find %d samples! Try to increase try_size.\n', d, size(samples,1));
                try_size = 100*batch_size;
            end
            if(isempty(idx))
                continue;
            end
            batch_num = batch_num + length(idx);
            samples = [samples; sample_try(idx,:)];
        end
        samples = samples(1:batch_size,:);
        sim_num = sim_num + batch_num;
        MCresults = isFailure(samples, threshold);
        if(isempty(find(MCresults,1)))
            R_min = R_tmp;
        else
            R_max = R_tmp;
        end
    end
    R_res = R_tmp;
end

