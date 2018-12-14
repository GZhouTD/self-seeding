clc;clear
close all

% foldname = 'elegant/scan/BC1MID/s14';
foldname = 'elegant/scan/BC1COL/dump_10M';
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
lambda = 2;
Wz = load(['elegant/wake/wiggler_wake_6per_',num2str(lambda),'um_filter.txt']);
tbin = Wz(:,1);
Wf = Wz(:,2);
t1 = tbin(tbin>0);
Wf = Wf(tbin>-max(t1));
tbin = tbin(tbin>-max(t1));
Wf = Wf*(zsep/6)*1;
dt = tbin(2) - tbin(1);

DL2R56 = 0e-6;
R56_chicane = -100e-6;
filter_coef = 3;

sddsfile = [foldname,'/col_3_s',num2str(ssnum),'/',plotname,'.out'];
plainfile = [foldname,'/col_3_s',num2str(ssnum),'/',plotname,'.6d'];
covert_command = ['sdds2plaindata ',sddsfile,' ',plainfile,' -col=x -col=xp -col=y -col=yp -col=t -col=p -sep=''  '''];
if ~exist(plainfile,'file')
    system(covert_command);
end
fid = fopen(plainfile);
np = fscanf(fid,'%d',1);
pdata = fscanf(fid,'%f',np*6);
fclose(fid);
x = pdata(1:6:np*6);
xp = pdata(2:6:np*6);
y = pdata(3:6:np*6);
yp = pdata(4:6:np*6);
t = pdata(5:6:np*6);
p = pdata(6:6:np*6);
t = t - mean(t);
slicepara = slice_beam(x,xp,y,yp,t,p,300,Q);

[X,Y,Z,~,~,I] = contour_plot_current(t,xp,300,300,0,Q);
figure
set(gcf,'position',[100,100,800,800])
subplot(2,1,1)
imagesc(X*1e15,Y,Z)
axis xy
colormap(jetvar)
xlabel('time (fs)')
ylabel('xp')
enhance_plot('times',16,2,8)
legend off
subplot(2,1,2)
plot(X*1e15,I)
xlabel('time (fs)')
ylabel('Current (kA)')
enhance_plot('times',16,2,8)
legend off

%%
%%%% beam matching parameters

tmatch1 = -4e-15;
tmatch2 = 3.2e-15;

x_core = x((t>tmatch1)&(t<tmatch2));
xp_core = xp((t>tmatch1)&(t<tmatch2));
y_core = y((t>tmatch1)&(t<tmatch2));
yp_core = yp((t>tmatch1)&(t<tmatch2));
p_core = p((t>tmatch1)&(t<tmatch2));
gamma0 = mean(p_core);

x = x - mean(x_core);
xp = xp - mean(xp_core);
y = y - mean(y_core);
yp = yp - mean(yp_core);


x_core = x_core - mean(x_core);
xp_core = xp_core - mean(xp_core);
y_core = y_core - mean(y_core);
yp_core = yp_core - mean(yp_core);

emitx0=sqrt(mean(x_core.^2).*mean(xp_core.^2)-mean(x_core.*xp_core).^2).*gamma0;
emity0=sqrt(mean(y_core.^2).*mean(yp_core.^2)-mean(y_core.*yp_core).^2).*gamma0;
betax0=mean(x_core.*x_core).*gamma0./emitx0;
betay0=mean(y_core.*y_core).*gamma0./emity0;
alphax0=-mean(x_core.*xp_core).*gamma0./emitx0;
alphay0=-mean(y_core.*yp_core).*gamma0./emity0;

clight = 299792458;
Lq = 0.3; % quad length
Ld = 3.6; % undulator section length
g = 12.64; % quad gradient
mc2 = gamma0*0.511e6;
kq = g*clight/mc2;
quadP = sqrt(kq)*Lq/2;

MF = [cos(quadP) sin(quadP)/sqrt(kq); -sqrt(kq)*sin(quadP) cos(quadP)];
MD = [cosh(quadP) sinh(quadP)/sqrt(kq); sqrt(kq)*sinh(quadP) cosh(quadP)];
ML = [1 Ld; 0 1];
1
A = MF*ML*MD*MD*ML*MF;
Cphi = A(1,1);
betaMAX = A(1,2)/sqrt(1-Cphi^2);

B = MD*ML*MF*MF*ML*MD;
Cphi = B(1,1);
betaMIN = B(1,2)/sqrt(1-Cphi^2);

mbetax = betaMAX;
mbetay = betaMIN;
malphax = 0;
malphay = 0;


%remove old correlation
xp=xp+alphax0.*x./betax0;
yp=yp+alphay0.*y./betay0;

%scale the beam
x = x.*sqrt(mbetax./betax0);
y = y.*sqrt(mbetay./betay0);
xp = xp.*sqrt(betax0./mbetax);
yp = yp.*sqrt(betay0./mbetay);

%add new correlation
xp = xp-malphax.*x./mbetax;
yp = yp-malphay.*y./mbetay;

x_core = x((t>tmatch1)&(t<tmatch2));
xp_core = xp((t>tmatch1)&(t<tmatch2));
y_core = y((t>tmatch1)&(t<tmatch2));
yp_core = yp((t>tmatch1)&(t<tmatch2));
p_core = p((t>tmatch1)&(t<tmatch2));

emitx1=sqrt(mean(x_core.^2).*mean(xp_core.^2)-mean(x_core.*xp_core).^2).*gamma0;
emity1=sqrt(mean(y_core.^2).*mean(yp_core.^2)-mean(y_core.*yp_core).^2).*gamma0;
betax1=mean(x_core.*x_core).*gamma0./emitx1;
betay1=mean(y_core.*y_core).*gamma0./emity1;
alphax1=-mean(x_core.*xp_core).*gamma0./emitx1;
alphay1=-mean(y_core.*yp_core).*gamma0./emity1;

[X,Y,Z,~,~,I] = contour_plot_current(t,xp,300,300,0,Q);
figure
set(gcf,'position',[100,100,800,800])
subplot(2,1,1)
imagesc(X*1e15,Y,Z)
axis xy
colormap(jetvar)
xlabel('time (fs)')
ylabel('xp')
enhance_plot('times',16,2,8)
legend off
subplot(2,1,2)
plot(X*1e15,I)
xlabel('time (fs)')
ylabel('Current (kA)')
enhance_plot('times',16,2,8)
legend off

%%

tcut1 = 5.5e-15;
pcut1 = 4986/0.511;

x((t<tcut1)&(p>pcut1)) = [];
xp((t<tcut1)&(p>pcut1)) = [];
y((t<tcut1)&(p>pcut1)) = [];
yp((t<tcut1)&(p>pcut1)) = [];
tbake = t;
pbake = p;
t((tbake<tcut1)&(pbake>pcut1)) = [];
p((tbake<tcut1)&(pbake>pcut1)) = [];
Q = Q/length(tbake)*length(t);
clear tbake pbake

[X,Y,Z,~,~,I] = contour_plot_current(t,xp,300,300,0,Q);
figure
set(gcf,'position',[100,100,800,800])
subplot(2,1,1)
imagesc(X*1e15,Y,Z)
axis xy
colormap(jetvar)
xlabel('time (fs)')
ylabel('xp')
enhance_plot('times',16,2,8)
legend off
subplot(2,1,2)
plot(X*1e15,I)
xlabel('time (fs)')
ylabel('Current (kA)')
enhance_plot('times',16,2,8)
legend off

delta = (p-mean(p))/mean(p);
t = t + DL2R56*delta/3e8;
[X,Y,Z,~,~,I] = contour_plot_current(t,p*0.511,300,300,0,Q);
figure(99)
set(gcf,'position',[100,100,800,800])
subplot(2,1,1)
imagesc(X*1e15,Y/1000,Z)
axis xy
colormap(jetvar)
xlabel('time (fs)')
ylabel('E (GeV)')
%     title(['no R56 effect'])
title([plotname])
%         set(gca,'xlim',[X(1)*1e15,X(end)*1e15]);
enhance_plot('times',16,2,8)
legend off
subplot(2,1,2)
% plot(slicepara.tbin/1e-15,slicepara.emitx)
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
figurefold = ['elegant/scan/BC1COL/dump_10M/figures/col_3_s',num2str(ssnum)];
if ~exist(figurefold,'dir')
    mkdir(figurefold);
end

for jj = 1:nsteps
    N = hist(t,tbin);
    I = Q/np*N/dt;
    I2 = I;
    I = filter1d(I,tbin,filter_coef*1e-15);
    I = abs(I);
    dE = conv(I*dt,Wf,'same');
    dE = dE';
    dp = interp1(tbin,dE/0.511,t);
    p = p + dp;
    delta = (p-mean(p))/mean(p);
    t = t + R56*delta/3e8+T566*delta.^2/3e8;
    t = t - mean(t);
    [X,Y,Z,~,~,I] = contour_plot_current(t,p*0.511,300,300,0,Q);
    figure(100)
    set(gcf,'position',[100,100,800,800])
    subplot(2,1,1)
    imagesc(X*1e15,Y/1000,Z)
    axis xy
    colormap(jetvar)
    xlabel('time (fs)')
    ylabel('E (GeV)')
    title(['no R56 effect'])
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
    if jj == nsteps
        figurename = [figurefold,'/first_chicane_',num2str(-R56_chicane/1e-6),'_filter_',num2str(filter_coef),'_und_',num2str(zsep*jj),'_',plotname,'.png'];
        print('-dpng','-r0',figurename)
    end
end
delta = (p-mean(p))/mean(p);
t = t + R56_chicane*delta/3e8;
[X,Y,Z,~,~,I] = contour_plot_current(t,p*0.511,300,300,0,Q);
figure(100)
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
enhance_plot('times',16,2,8)
legend off
figurename = [figurefold,'/chicane_',num2str(-R56_chicane/1e-6),'_filter_',num2str(filter_coef),'_und_',num2str(zsep*jj),'_',plotname,'.png'];
print('-dpng','-r0',figurename)
%     pause
for jj = 1:nsteps2
    N = hist(t,tbin);
    I = Q/np*N/dt;
    I = filter1d(I,tbin,filter_coef*1e-15);
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
    prange = [4.94e3,5.01e3];
    %         prange = [3.4e3,3.5e3];
    [X,Y,Z,~,~,I] = contour_plot_current(t,p*0.511,300,300,0,Q,trange,prange);
    figure(100)
    set(gcf,'position',[100,100,800,800])
    subplot(2,1,1)
    imagesc(X*1e15,Y/1000,Z)
    axis xy
    colormap(jetvar)
    xlabel('time (fs)')
    ylabel('E (GeV)')
    title(['no R56 effect'])
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
    enhance_plot('times',16,2,8)
    legend off
    %         figurename = [figurefold,'/no_R56.png'];
    if jj == nsteps2
        figurename = [figurefold,'/second_chicane_',num2str(-R56_chicane/1e-6),'_filter_',num2str(filter_coef),'_und_',num2str(zsep*jj),'_',plotname,'.png'];
        print('-dpng','-r0',figurename)
    end
end


[X,Y,Z,~,~,I] = contour_plot_current(t,xp,300,300,0,Q);
figure
set(gcf,'position',[100,100,800,800])
subplot(2,1,1)
imagesc(X*1e15,Y,Z)
axis xy
colormap(jetvar)
xlabel('time (fs)')
ylabel('xp')
enhance_plot('times',16,2,8)
legend off
subplot(2,1,2)
plot(X*1e15,I)
xlabel('time (fs)')
ylabel('Current (kA)')
enhance_plot('times',16,2,8)
legend off


pdata_new = [x,xp,y,yp,t,p];
dumpfile = [foldname,'/col_3_s',num2str(ssnum),'/',plotname,'_chicane_',num2str(-R56_chicane/1e-6),'_filter_',num2str(filter_coef),'_match2_.wig'];
if ~exist(dumpfile,'file')
    fid3 = fopen(dumpfile,'w');
    fprintf(fid3, '%12.10E\t %12.10E\t %12.10E\t %12.10E\t %12.10E\t %12.10E\t\n',pdata_new');
    fclose(fid3);
    fclose all;
end