clear
fid = fopen('MCres.txt', 'r');
MCpfail = [];
MCfom=[];
sample_n = 1000;

while(1)
    line = fgetl(fid);
    if(~ischar(line))
        break;
    end
    if(strfind(line,'failed'))
        %line = fgetl(fid);
        id_s=strfind(line,'MC failure rate = ')+18;
        id_e=strfind(line,'MC FOM = ') + 9;
        line1 = line(id_s:id_e-12);
        MCresults = sscanf(line1, '%g');
        MCpfail = [MCpfail, MCresults];
        line2 = line(id_e:end);
        MCresults = sscanf(line2, '%g');
        MCfom = [MCfom, MCresults];
        sample_n = [sample_n, sample_n(end)+1000];
    end
end
fclose(fid);
figure(1)
semilogx(sample_n(2:end),MCpfail, '-*'); title('MCPFail');
xlabel('number of samples');
figure(2)
semilogx(sample_n(2:end),MCfom, '-*'); title('MCfom');
xlabel('number of samples');
