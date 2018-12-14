clc;clear
close all

foldname = 'elegant/scan/BC1COL/col_3';
ssnum = 59;
zsep = 1;
nwigg = 10;
nwigg2 = 24;

x_max_col = 3;

charge = 189.4225e-12; % x_max_col = 3;

Q = charge;
plotname = 'XLBEG';
lambda = 3;
Wz = load(['elegant/wake/wiggler_wake_6per_',num2str(lambda),'um_filter.txt']);
tbin = Wz(:,1);
Wf = Wz(:,2);
t1 = tbin(tbin>0);
Wf = Wf(tbin>-max(t1));
tbin = tbin(tbin>-max(t1));
Wf = Wf*(zsep/6)*1;
dt = tbin(2) - tbin(1);

DL2R56 = 0e-6;
R56_chicane = -250e-6;

for ii = 1:length(ssnum)
    sddsfile = [foldname,'/s',num2str(ssnum(ii)),'/',plotname,'.out'];
    plainfile = [foldname,'/s',num2str(ssnum(ii)),'/',plotname,'.pla'];
    covert_command = ['sdds2plaindata ',sddsfile,' ',plainfile,' -col=t -col=p -sep=''  '''];
    if ~exist(plainfile,'file')
        system(covert_command);
    end
    fid = fopen(plainfile);
    np = fscanf(fid,'%d',1);
    pdata = fscanf(fid,'%f',np*2);
    fclose(fid);
    t = pdata(1:2:np*2);
    p = pdata(2:2:np*2);
    t = t - mean(t);
    delta = (p-mean(p))/mean(p);
    t = t + DL2R56*delta/3e8;
    [X,Y,Z,~,~,I] = contour_plot_current(t,p*0.511,300,300,0,Q);
    figure(2)
    set(gcf,'position',[100,100,800,800])
    subplot(2,1,1)
    imagesc(X*1e15,Y/1000,Z)
    axis xy
    colormap(jetvar)
    xlabel('time (fs)')
    ylabel('E (GeV)')
    title([plotname])
    enhance_plot('times',16,2,8)
    legend off
    subplot(2,1,2)
    plot(X*1e15,I)
    xlabel('time (fs)')
    ylabel('Current (kA)')
    set(gca,'xlim',[X(1)*1e15,X(end)*1e15]);
    set(gca,'ylim',[0, max(I)*1.2]);
    enhance_plot('times',16,2,8)
    legend off
    nsteps = nwigg/zsep;
    nsteps2 = nwigg2/zsep;
    R56 = -2*lambda*1e-6*zsep;
    T566 = -1.5*R56;
    figurefold = ['elegant/scan/figures/BC1COL/col_3/s',num2str(ssnum(ii))];
    if ~exist(figurefold,'dir')
        mkdir(figurefold);
    end
    
    for jj = 1:nsteps
        N = hist(t,tbin);
        I = Q/np*N/dt;
        I2 = I;
        I = filter1d(I,tbin,4e-15);
        I = abs(I);
        dE = conv(I*dt,Wf,'same');
        dE = dE';
        dp = interp1(tbin,dE/0.511,t);
        p = p + dp;
        delta = (p-mean(p))/mean(p);
        t = t + R56*delta/3e8+T566*delta.^2/3e8;
        t = t - mean(t);
        [X,Y,Z,~,~,I] = contour_plot_current(t,p*0.511,300,300,0,Q);
        figure(1)
        set(gcf,'position',[100,100,800,800])
        subplot(2,1,1)
        imagesc(X*1e15,Y/1000,Z)
        axis xy
        colormap(jetvar)
        xlabel('time (fs)')
        ylabel('E (GeV)')
        title([plotname,': ZSEP=',num2str(zsep),'  Und. #=',num2str(zsep*jj)])
        enhance_plot('times',16,2,8)
        legend off
        subplot(2,1,2)
        plot(X*1e15,I)
        xlabel('time (fs)')
        ylabel('Current (kA)')
        set(gca,'xlim',[X(1)*1e15,X(end)*1e15]);
        set(gca,'ylim',[0, max(I)*1.2]);
        enhance_plot('times',16,2,8)
        legend off
    end
    [tt,pavg,prms] = sliceEnergy(t,p,300);
    figure
    plot(tt,pavg*0.511)
    xlabel('t (fs)')
    ylabel('p_{avg}')
    set(gca,'ylim',[4980,5020])
    enhance_plot('times',16,2,8)
    yyaxis right
    plot(tt,prms*0.511)
    set(gca,'ylim',[0,5])
    ylabel('p_{rms}')
    legend('mean energy','rms energy')
    enhance_plot('times',16,2,8)
    
%     pause
%     delta = (p-mean(p))/mean(p);
%     t = t + R56_chicane*delta/3e8;
%     [X,Y,Z,~,~,I] = contour_plot_current(t,p*0.511,300,300,0,Q);
%     figure(1)
%     set(gcf,'position',[100,100,800,800])
%     subplot(2,1,1)
%     imagesc(X*1e15,Y/1000,Z)
%     axis xy
%     colormap(jetvar)
%     xlabel('time (fs)')
%     ylabel('E (GeV)')
%     title([plotname,': chicane'])
%     enhance_plot('times',16,2,8)
%     legend off
%     subplot(2,1,2)
%     plot(X*1e15,I)
%     xlabel('time (fs)')
%     ylabel('Current (kA)')
%     set(gca,'xlim',[X(1)*1e15,X(end)*1e15]);
%     set(gca,'ylim',[0, max(I)*1.2]);
%     enhance_plot('times',16,2,8)
%     legend off
    
end