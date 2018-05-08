function DBIres = DBI(data, idx, C)

[cluster_num,~] = size(C);
CP = zeros(cluster_num,1);
for i = 1:cluster_num
    data_k = data(idx==i,:);
    C_k = C(i,:);
    CP(i) = 1/size(data_k,1)* ...
        sum(Cluster_norm(C_k, data_k));
end

DBIres = 0;
for i = 1:cluster_num
    tmp = zeros(cluster_num,1);
    for j = 1:cluster_num
        if(i==j)
            continue;
        end
        tmp(j) = (CP(i)+ CP(j))/(Cluster_norm(C(i,:),C(j,:)));
    end
    DBIres = DBIres + max(tmp);
end
DBIres = 1/cluster_num * DBIres;

if(DBIres==0)
    disp('error');
end

end

