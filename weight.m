function res = weight(samples,C,alpha,beta)

cluster_num=size(C,1);
[N,~] = size(samples);

fx = mvnpdf(samples);
res = zeros(N,1);
for i = 1:cluster_num
    shifted_samples = bsxfun(@minus, samples, C(i,:));
    res = res + beta(i) * (mvnpdf(shifted_samples)./fx);
end

res = alpha+(1-alpha)*res;
res = 1./res;

end

