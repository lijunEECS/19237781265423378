function res = isFailure(data)

data_norm = sqrt(sum(data.^2,2));

complex_data = data(:,1) + 1i*data(:,2);

data_angle = angle(complex_data);

res1 = bsxfun(@and, data_norm>4.8, data_angle >= 2*pi/3);
res1 = bsxfun(@and, res1, data_angle <= 3*pi/4);

res2 = bsxfun(@and, data_norm>3.9, data_angle >= -2*pi/3);
res2 = bsxfun(@and, res2, data_angle <= -pi/2);

res3 = bsxfun(@and, data_norm>4.9, data_angle >= pi/3);
res3 = bsxfun(@and, res3, data_angle <= pi/2);

res4 = bsxfun(@and, data_norm>3.5, data_angle >= -pi/5);
res4 = bsxfun(@and, res4, data_angle <= 0);

res = bsxfun(@or, res1, res2);
res = bsxfun(@or, res, res3);
res = bsxfun(@or, res, res4);

end

