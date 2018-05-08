function w = sampleWeight(samples, R_min)

% data = samples(:,1:end-1);
% [~,D] = size(data);
% 
% data1 = data(:,1:end/2);
% data2 = data(:,end/2+1:end);
% p1 = mvnpdf(data1,zeros(1,D/2),ones(1,D/2));
% p2 = mvnpdf(data2,zeros(1,D/2),ones(1,D/2));
% 
% w1 = 1/sum(p1)*p1;
% w2 = 1/sum(p2)*p2;
% 
% w = w1.*w2;
% w = 1/sum(w)*w;

R = samples(:,end);

w = exp(-(1/R_min)*R);

w = w/sum(w);

end

