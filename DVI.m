function DVIres = DVI(data, idx)

cluster_num = max(idx);
max_dis_inC = zeros(cluster_num,1);

for i = 1:cluster_num
    data_k = data(idx==i,:);
    N = size(data_k,1);
    max_dis = 0;
    for j = 1:N
        for k = 1:N
            if(j==k)
                continue;
            end
            tmp_dis = Cluster_norm(data_k(j,:),data_k(k,:));
            if(tmp_dis>max_dis)
                max_dis = tmp_dis;
            end
        end
    end
    max_dis_inC(i) = max_dis;
end

min_dis_betC = 5e150*ones(cluster_num,cluster_num);
for i = 1:cluster_num
    for j = 1:cluster_num
        if(i==j)
            continue;
        end
        data_i = data(idx==i,:);
        data_j = data(idx==j,:);
        N_i = size(data_i,1);
        N_j = size(data_j,1);
        min_dis = 5e150;
        for p=1:N_i
            for q=1:N_j
                tmp_dis = Cluster_norm(data_i(p,:),data_j(q,:));
                if(tmp_dis<min_dis)
                    min_dis = tmp_dis;
                end
            end
        end
        min_dis_betC(i,j) = min_dis;
    end
end

DVIres=min(min(min_dis_betC))/max(max_dis_inC);

end