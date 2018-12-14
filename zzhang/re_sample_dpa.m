clc;clear
close all
filename = 'genesis/XLBEG2_ssnum_51_chicane_100_filter_3/s1/genesis.out';
[data,info] = readOutput(filename);
dpa = gen_load_dpa(filename);


%%
N = 1e6;
nslice = info.nt;
cur = data.cur;
cur_max = max(cur);
cur_sum = sum(cur);
np = 16384;
sq = [];
for is = 1:length(cur)
    npart = round(N*cur(is)/cur_sum);
    dn = round(np/npart);
    sq = [sq,((is-1)*np+1):dn:(is*np)];
end

pdata = dpa(sq,:);

t = pdata(:,5);
p = pdata(:,6);

contour_plot(t/3e8,p,300,300,1,69.3669e-12);

% fid3 = fopen(distfileout,'w');
% fprintf(fid3, '%12s\n', '# Double-Taper Attosection XFEL');
% fprintf(fid3, '%12s\n', '? version = 1.0');
% fprintf(fid3, '%11s  %e\n', '? charge = ',Qtot/length(dataSl(:,1))*length(newdata(:,1)));
% fprintf(fid3, '%9s  %d\n', '? size = ',length(newdata(:,1)));
% fprintf(fid3, '%12s\n', '? COLUMNS X XPRIME Y YPRIME Z GAMMA');
% 
% fprintf(fid3, '%12.10E\t %12.10E\t %12.10E\t %12.10E\t %12.10E\t %12.10E\t\n',newdata');
% fclose(fid3);
% 
% fclose all;