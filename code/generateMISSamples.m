function [samples, w] = generateMISSamples(C, alpha, beta, batch_size)

[cluster_num,D] = size(C);

samples = mvnrnd(zeros(1,D),ones(1,D),batch_size);

alpha_rnd = rand([batch_size,1]);
beta_rnd = rand([batch_size,1]);
acc_beta = cumsum(beta);
C_id = sum(bsxfun(@gt,repmat(beta_rnd',[cluster_num,1]),acc_beta))+1;

for i=1:batch_size
    if alpha_rnd<alpha
        continue;
    else
        samples(i,:) = samples(i,:) + C(C_id(i),:);
    end
end

w = weight(samples, C, alpha, beta);

end

