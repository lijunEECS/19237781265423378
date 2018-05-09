stop_fom = 0.1;

MCpfail = [];
MCfom = [];
sample_n = 0;
D = 360;

sample_unit = 100;
iter = 1;
threshold = 135;
MCtotal_error_counter = 0;

while(1)
    samples = normrnd(0,1,sample_unit, D);
    MCresults = isFailure(samples, threshold);
    MCerror_counter = sum(MCresults);
    MCtotal_error_counter = MCtotal_error_counter + MCerror_counter;
    sample_n = [sample_n, sample_n(end) + sample_unit];
    MCpfail = [MCpfail, MCtotal_error_counter/sample_n(end)];
    MCfom = [MCfom, std(MCpfail)/mean(MCpfail)];
    str = sprintf('%d out of %d samples failed(%d/%d), MC failure rate = %e, MC FOM = %e', MCtotal_error_counter, sample_n(end), MCerror_counter, sample_unit, MCpfail(end), MCfom(end));
    disp(str);
    if(MCfom(end)<stop_fom && iter>5)
        break;
    end
    iter = iter+1;
end