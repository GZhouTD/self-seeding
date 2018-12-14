clc;clear
close all

foldname = 'genesis/dump_1_XLBEG_ssnum_51_chicane_100_filter_3';
clight = 2.998e8;
R56 = 10e-15*clight;
ssnum = 51:100;

for ii = 1:length(ssnum)
    ii
    filename = [foldname,'/s',num2str(ssnum(ii)),'/genesis.out'];
    [data,info] = readOutput(filename);
    dpa = gen_load_dpa(filename);
    
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
    
    t_cut = pdata(:,5);
    p_cut = pdata(:,6);
    x_cut = pdata(:,1);
    xp_cut = pdata(:,2);
    y_cut = pdata(:,3);
    yp_cut = pdata(:,4);
    delta = (p_cut - mean(p_cut))/mean(p_cut);
    t_cut = t_cut + R56 * delta;
    t_cut = t_cut - min(t_cut);
    t_cut = t_cut/clight;
    Q = sum(data.cur.*info.lambda/clight);
    contour_plot(t_cut,p_cut,300,300,1,Q);
    
    pdata_new = [x_cut,xp_cut,y_cut,yp_cut,t_cut*clight,p_cut];
    
    distfile = [foldname,'/s',num2str(ssnum(ii)),'/beam2.dist'];
    
    fid3 = fopen(distfile,'w');
    fprintf(fid3, '%12s\n', '# Double-Taper Attosection XFEL');
    fprintf(fid3, '%12s\n', '? version = 1.0');
    fprintf(fid3, '%11s  %e\n', '? charge = ',Q);
    fprintf(fid3, '%9s  %d\n', '? size = ',length(t_cut));
    fprintf(fid3, '%12s\n', '? COLUMNS X XPRIME Y YPRIME Z GAMMA');
    
    fprintf(fid3, '%12.10E\t %12.10E\t %12.10E\t %12.10E\t %12.10E\t %12.10E\t\n',pdata_new');
    fclose(fid3);
    
    fclose all;
end