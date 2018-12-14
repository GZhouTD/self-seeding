clc;clear
close all

lambda = 3;
zsep = 1;
nsteps = 7;
Wz = load(['elegant/wake/wiggler_wake_6per_',num2str(lambda),'um_filter.txt']);
tbin = Wz(:,1);
Wf = Wz(:,2);
t1 = tbin(tbin>0);
Wf = Wf(tbin>-max(t1));
tbin = tbin(tbin>-max(t1));
Wf = Wf*(zsep/6)*1;
dt = tbin(2) - tbin(1);
R56 = -2*lambda*1e-6*zsep;
T566 = -1.5*R56;
R56_xleap = -400e-6;

np = 1e6;
sig_t = 3.5e-15/2.35;
gamma0 = 4.5e3/0.511;
sig_gamma = 1/0.511/6*3;
t1 = sig_t*randn(np,1);
gamma1 = 2*sig_gamma*randn(np,1)+gamma0;

np2 = 2.3e6;
Lt = 30e-15/6*3.5;
t2 = -Lt/2+Lt*rand(np2,1);
gamma2 = sig_gamma*randn(np2,1)+gamma0;

p = [gamma1;gamma2];
t = [t1;t2];
np = length(t);
Q = 85e-12/6*3.5;
[X,Y,Z,~,~,I] = contour_plot_current(t,p*0.511,300,300,0,Q);

pp = [];
sumI = sum(I);
t = sort(t,'ascend');
for ii = 1:length(I)
    nn = round(I(ii)./sumI*np);
    pp_slice = sig_gamma*I(ii)*randn(nn,1)+gamma0;
    pp = [pp;pp_slice];
end
p = pp;
t = t(1:length(p));
[X,Y,Z,~,~,I] = contour_plot_current(t,p*0.511,300,300,0,Q);

figure(1)
% set(gcf,'position',[100,100,800,800])
subplot(2,1,1)
imagesc(X*1e15,Y/1000,Z)
axis xy
colormap(jetvar)
xlabel('time (fs)')
ylabel('E (GeV)')
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

%%

for jj = 1:nsteps
    N = hist(t,tbin);
    I = Q/np*N/dt;
    I2 = I;
    I = filter1d(I,tbin,0.25e-15);
    I = abs(I);
    dE = conv(I*dt,Wf,'same');
    dE = dE';
    dp = interp1(tbin,dE/0.511,t);
    p = p + dp;
    delta = (p-mean(p))/mean(p);
    t = t + R56*delta/3e8+T566*delta.^2/3e8;
    t = t - mean(t);
    
end

[X,Y,Z,~,~,I] = contour_plot_current(t,p*0.511,300,300,0,Q);
figure(2)
subplot(2,1,1)
imagesc(X*1e15,Y/1000,Z)
axis xy
colormap(jetvar)
xlabel('time (fs)')
ylabel('E (GeV)')
title([': ZSEP=',num2str(zsep),'  Und. #=',num2str(zsep*jj)])
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


delta = (p-mean(p))/mean(p);
t = t + R56_xleap*delta/3e8;
t = t - mean(t);
[X,Y,Z,~,~,I] = contour_plot_current(t,p*0.511,300,300,0,Q);
figure(3)
subplot(2,1,1)
imagesc(X*1e15,Y/1000,Z)
axis xy
colormap(jetvar)
xlabel('time (fs)')
ylabel('E (GeV)')
title(['after xleap chicane'])
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