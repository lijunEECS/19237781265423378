function idx = reindex(idx)

N = max(idx);

idx_cnt = zeros(N,1);
for i=1:N
    idx_cnt(i) = sum(idx==i);
end

for i=2:N
    last_nz_id = i-1;
    while(idx_cnt(last_nz_id)==0)
        last_nz_id = last_nz_id-1;
        if(last_nz_id==0)
            break;
        end
    end
    idx(idx==i) = last_nz_id+1;
    if(last_nz_id ~= i-1)
        idx_cnt(last_nz_id+1) = idx_cnt(i);
        idx_cnt(i)=0;
    end
end

end

