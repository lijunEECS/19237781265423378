function [idx, C, k] = Wkmeans(samples, k, w)

sample_num = size(samples, 1);

idx = zeros(sample_num, 1);
last_idx = ones(sample_num, 1);

C_idx = randperm(sample_num, k);
C = samples(C_idx,:);

while(idx~=last_idx)
    last_idx = idx;
    [~,idx] = max(samples*C', [], 2);
    
    idx_cnt = zeros(k,1);
    for i=1:k
        idx_cnt(i) = sum(idx==i);
    end
    
    rm_idx = find(idx_cnt==0,1);
    
    if(~isempty(rm_idx))
        N = length(rm_idx);
        if(k<=N)
            disp('error');
        end
        k = k - N;
        for i=1:N
            C(rm_idx(i),:)=[];
        end
        idx=reindex(idx);
    end
    
    for i=1:k
        samples_k = samples(idx==i,:);
        w_k = w(idx==i);
        tmp = bsxfun(@times,samples_k,w_k);
        C(i,:) = sum(tmp,1);
        C(i,:) = 1/norm(C(i,:))*C(i,:);
    end
end

if(k==0)
    disp('error');
end

end

