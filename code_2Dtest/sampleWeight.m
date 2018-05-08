function w = sampleWeight(samples, R_min)

R=samples(:,end);

w = exp(-(1/R_min)*R);

w = w/sum(w);

end

