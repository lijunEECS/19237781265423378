function d = Cluster_norm(C, samples)

C_norm = norm(C);
samples_num = size(samples,1);
d = zeros(samples_num,1);
for i = 1:samples_num
    d(i) = 1 - (1/(norm(samples(i,:))*C_norm))*(samples(i,:)*C');
end

end

