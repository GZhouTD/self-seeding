clc;clear
close all

% foldname = 'elegant/scan/BC1MID/s14';
foldname = 'elegant/scan/BC1COL/col_3';
ssnum = 51;
zsep = 1;
nwigg = 18;
nwigg2 = 24;

x_max_col = 3;

% charge = 249.2432e-12;  % x_max_col = 5;
% charge = 231.9349e-12; % x_max_col = 4;
% charge =  217.5695e-12; % x_max_col = 3.6;
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
R56_chicane = -200e-6;

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
%     slicepara = 
    figure(2)
    set(gcf,'position',[100,100,800,800])
    subplot(2,1,1)
    imagesc(X*1e15,Y/1000,Z)
    axis xy
    colormap(jetvar)
    xlabel('time (fs)')
    ylabel('E (GeV)')
    pause
%     title(['no R56 effect'])
    title([plotname])
    %         set(gca,'xlim',[X(1)*1e15,X(end)*1e15]);
    enhance_plot('times',16,2,8)
    legend off
    subplot(2,1,2)
%     hold on
%     yyaxis right
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
%         figure(101)
%         plot(I2)
%         hold on
%         plot(I);
%         hold off
        dE = conv(I*dt,Wf,'same');
%         dE = smooth1d(dE,4);
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
%         title(['no R56 effect'])
        title([plotname,': ZSEP=',num2str(zsep),'  Und. #=',num2str(zsep*jj)])
        %         set(gca,'xlim',[X(1)*1e15,X(end)*1e15]);
        enhance_plot('times',16,2,8)
        legend off
        subplot(2,1,2)
        plot(X*1e15,I)
        xlabel('time (fs)')
        ylabel('Current (kA)')
        set(gca,'xlim',[X(1)*1e15,X(end)*1e15]);
        set(gca,'ylim',[0, max(I)*1.2]);
%         title([plotname,': ZSEP=',num2str(zsep),'  Und. #=',num2str(zsep*jj)])
        enhance_plot('times',16,2,8)
        legend off
        %         figurename = [figurefold,'/no_R56.png'];
%         figurename = [figurefold,'/zsep_',num2str(zsep),'_und_',num2str(zsep*jj),'_',plotname,'_T566.png'];
        %         print('-dpng','-r0',figurename)
    end
    pause
    delta = (p-mean(p))/mean(p);
    t = t + R56_chicane*delta/3e8;
    [X,Y,Z,~,~,I] = contour_plot_current(t,p*0.511,300,300,0,Q);
    figure(1)
    set(gcf,'position',[100,100,800,800])
    subplot(2,1,1)
    imagesc(X*1e15,Y/1000,Z)
    axis xy
    colormap(jetvar)
    xlabel('time (fs)')
    ylabel('E (GeV)')
%     title(['no R56 effect'])
    title([plotname,': chicane'])
    %         set(gca,'xlim',[X(1)*1e15,X(end)*1e15]);
    enhance_plot('times',16,2,8)
    legend off
    subplot(2,1,2)
    plot(X*1e15,I)
    xlabel('time (fs)')
    ylabel('Current (kA)')
    set(gca,'xlim',[X(1)*1e15,X(end)*1e15]);
    set(gca,'ylim',[0, max(I)*1.2]);
%     title([plotname,': chicane'])
    enhance_plot('times',16,2,8)
    legend off
    pause
    for jj = 1:nsteps2
        N = hist(t,tbin);
        I = Q/np*N/dt;
        I = filter1d(I,tbin,4e-15);
        I = abs(I);
        dE = conv(I*dt,Wf,'same');
%         dE = smooth1d(dE,8);
        dE = dE';
        dp = interp1(tbin,dE/0.511,t);
        p = p + dp;
        delta = (p-mean(p))/mean(p);
        t = t + R56*delta/3e8+T566*delta.^2/3e8;
        t = t - mean(t);
        trange = [-40e-15,20e-15];
%         prange = [4.97e3,5.03e3];
        prange = [4.92e3,5.01e3];
%         prange = [3.4e3,3.5e3];
        [X,Y,Z,~,~,I] = contour_plot_current(t,p*0.511,300,300,0,Q,trange,prange);
        figure(1)
        set(gcf,'position',[100,100,800,800])
        subplot(2,1,1)
        imagesc(X*1e15,Y/1000,Z)
        axis xy
        colormap(jetvar)
        xlabel('time (fs)')
        ylabel('E (GeV)')
%         title(['no R56 effect'])
        title([plotname,': ZSEP=',num2str(zsep),'  Und. #=',num2str(zsep*jj)])
        %         set(gca,'xlim',[X(1)*1e15,X(end)*1e15]);
        enhance_plot('times',16,2,8)
        legend off
        subplot(2,1,2)
        plot(X*1e15,I)
        xlabel('time (fs)')
        ylabel('Current (kA)')
        set(gca,'xlim',[X(1)*1e15,X(end)*1e15]);
        set(gca,'ylim',[0, max(I)*1.2]);
%         title([plotname,': ZSEP=',num2str(zsep),'  Und. #=',num2str(zsep*jj)])
        enhance_plot('times',16,2,8)
        legend off
        %         figurename = [figurefold,'/no_R56.png'];
%         figurename = [figurefold,'/zsep_',num2str(zsep),'_und_',num2str(zsep*jj),'_',plotname,'_T566.png'];
        %         print('-dpng','-r0',figurename)
    end
    
end