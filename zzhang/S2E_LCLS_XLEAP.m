clc;clear
close all
% fid = fopen('elegant/s0_10M/pdata');
% np = 11976490;
fid = fopen('elegant/scan/BC1COL/col_3/s51/XLBEG.pla');
np = fscanf(fid,'%d',1);

% fid = fopen('elegant/s3_1M/pdata');
% np = fscanf(fid,'%d',1);

pdata = fscanf(fid,'%f',np*6);
fclose(fid);
ngroup = floor(np/6);
sq = (1:ngroup)*6-6;
x = pdata(sq+1);
xp = pdata(sq+2);
y = pdata(sq+3);
yp = pdata(sq+4);
t = pdata(sq+5);
p = pdata(sq+6);
Np = length(x);

%%
sliceQ = 0;
dumpQ = 0;
t = t - mean(t);

charge = 189.4225e-12; % x_max_col = 3;

Q = charge;


% p_cut_min1 = 6792;
% 
% x(p>p_cut_min1) = [];
% xp(p>p_cut_min1) = [];
% y(p>p_cut_min1) = [];
% yp(p>p_cut_min1) = [];
% t(p>p_cut_min1) = [];
% p(p>p_cut_min1) = [];


[X,Y,Z,dx,dy] = contour_plot(t,p,300,300,0,Q);

figure
imagesc(X*1e15,Y,Z)
axis xy
xlabel('time (fs)')
ylabel('\gamma')
colormap(jetvar)
set(gca,'xlim',[-30 30])
enhance_plot('times',16,2,8)
legend off

Wz = load('elegant/wake/wiggler_wake_6per_4um.txt');
Wf_factor = 4;

tbin = Wz(:,1);
Wf = Wz(:,2)*Wf_factor;
t1 = tbin(tbin>0);
Wf = Wf(tbin>-max(t1));
tbin = tbin(tbin>-max(t1));

dt = tbin(2) - tbin(1);

N = hist(t,tbin);
% Q = 250e-12;
I = Q/np*N/dt;

slicepara = slice_beam(x,xp,y,y,t,p,300,Q);

figure
plot(tbin/1e-15,I)



figure
plot(Wz(:,1)/1e-15,Wz(:,2))

dE = conv(I*dt,Wf,'same');
dE = smooth1d(dE);
dE = smooth1d(dE);
dE = dE';

dp = interp1(tbin,dE/0.511,t);
p = p + dp;
[X,Y,Z,dx,dy] = contour_plot(t,p,300,300,0,Q);

I = sum(Z,1)/dx*Q/Np;

figure
imagesc(X*1e15,Y,Z)
axis xy
xlabel('time (fs)')
ylabel('\gamma')
colormap(jetvar)
set(gca,'xlim',[-30 30])
enhance_plot('times',16,2,8)
legend off

figure
plot(X*1e15,I/1000)
xlabel('time (fs)')
ylabel('Current (kA)')
set(gca,'xlim',[-30 30])
enhance_plot('times',16,2,8)
legend off
xlamd = 3e-2;
gamma0 = mean(p);
Krms = 2.5;
xlamds = (1+Krms^2)/2/gamma0^2*xlamd;


%%

%%%% beam matching parameters
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



t_cut_min2 = -5.5e-15;
t_cut_max2 = 9.5e-15;
p_cut_min2 = 6722;
p_cut_max2 = 6880;

x_ecut2 = x((p>p_cut_min2)&(p<p_cut_max2)&(t>t_cut_min2)&(t<t_cut_max2));
xp_ecut2 = xp((p>p_cut_min2)&(p<p_cut_max2)&(t>t_cut_min2)&(t<t_cut_max2));
y_ecut2 = y((p>p_cut_min2)&(p<p_cut_max2)&(t>t_cut_min2)&(t<t_cut_max2));
yp_ecut2 = yp((p>p_cut_min2)&(p<p_cut_max2)&(t>t_cut_min2)&(t<t_cut_max2));
t_ecut2 = t((p>p_cut_min2)&(p<p_cut_max2)&(t>t_cut_min2)&(t<t_cut_max2));
p_ecut2 = p((p>p_cut_min2)&(p<p_cut_max2)&(t>t_cut_min2)&(t<t_cut_max2));

t_ecut2 = -t_ecut2;
t_ecut2 = t_ecut2 - min(t_ecut2);

contour_plot(t_ecut2,p_ecut2,300,300,1);

t_cut_min3 = 11.6e-15;
p_cut_min3 = 6729;

x_ecut2((p_ecut2<p_cut_min3)&(t_ecut2<t_cut_min3)) = [];
xp_ecut2((p_ecut2<p_cut_min3)&(t_ecut2<t_cut_min3)) = [];
y_ecut2((p_ecut2<p_cut_min3)&(t_ecut2<t_cut_min3)) = [];
yp_ecut2((p_ecut2<p_cut_min3)&(t_ecut2<t_cut_min3)) = [];
t_ecut2_cp = t_ecut2;
p_ecut2_cp = p_ecut2;
t_ecut2((p_ecut2_cp<p_cut_min3)&(t_ecut2_cp<t_cut_min3)) = [];
p_ecut2((p_ecut2_cp<p_cut_min3)&(t_ecut2_cp<t_cut_min3)) = [];

t_cut_min3 = 10.93e-15;
p_cut_min3 = 6755;

x_ecut2((p_ecut2<p_cut_min3)&(t_ecut2<t_cut_min3)) = [];
xp_ecut2((p_ecut2<p_cut_min3)&(t_ecut2<t_cut_min3)) = [];
y_ecut2((p_ecut2<p_cut_min3)&(t_ecut2<t_cut_min3)) = [];
yp_ecut2((p_ecut2<p_cut_min3)&(t_ecut2<t_cut_min3)) = [];
t_ecut2_cp = t_ecut2;
p_ecut2_cp = p_ecut2;
t_ecut2((p_ecut2_cp<p_cut_min3)&(t_ecut2_cp<t_cut_min3)) = [];
p_ecut2((p_ecut2_cp<p_cut_min3)&(t_ecut2_cp<t_cut_min3)) = [];

t_cut_min3 = 5e-15;
t_cut_max3 = 9.6e-15;
p_cut_min3 = 6771;

x_ecut2((p_ecut2<p_cut_min3)&(t_ecut2<t_cut_max3)&(t_ecut2>t_cut_min3)) = [];
xp_ecut2((p_ecut2<p_cut_min3)&(t_ecut2<t_cut_max3)&(t_ecut2>t_cut_min3)) = [];
y_ecut2((p_ecut2<p_cut_min3)&(t_ecut2<t_cut_max3)&(t_ecut2>t_cut_min3)) = [];
yp_ecut2((p_ecut2<p_cut_min3)&(t_ecut2<t_cut_max3)&(t_ecut2>t_cut_min3)) = [];
t_ecut2_cp = t_ecut2;
p_ecut2_cp = p_ecut2;
t_ecut2((p_ecut2_cp<p_cut_min3)&(t_ecut2_cp<t_cut_max3)&(t_ecut2_cp>t_cut_min3)) = [];
p_ecut2((p_ecut2_cp<p_cut_min3)&(t_ecut2_cp<t_cut_max3)&(t_ecut2_cp>t_cut_min3)) = [];


% t_cut_min3 = 1.23e-15;
% p_cut_min3 = 6860;
% 
% x_ecut2((p_ecut2>p_cut_min3)&(t_ecut2>t_cut_min3)) = [];
% xp_ecut2((p_ecut2>p_cut_min3)&(t_ecut2>t_cut_min3)) = [];
% y_ecut2((p_ecut2>p_cut_min3)&(t_ecut2>t_cut_min3)) = [];
% yp_ecut2((p_ecut2>p_cut_min3)&(t_ecut2>t_cut_min3)) = [];
% t_ecut2_cp = t_ecut2;
% p_ecut2_cp = p_ecut2;
% t_ecut2((p_ecut2_cp>p_cut_min3)&(t_ecut2_cp>t_cut_min3)) = [];
% p_ecut2((p_ecut2_cp>p_cut_min3)&(t_ecut2_cp>t_cut_min3)) = [];


z = clight*t_ecut2;
[X,Y,Z,dx,dy] = contour_plot(t_ecut2,p_ecut2,300,300,0);
I = sum(Z,1)/dx*Q/Np;


figure
yyaxis left
imagesc(X/1e-15,Y,Z)
axis xy
xlabel('time (fs)')
ylabel('\gamma')
colormap(jetvar)
enhance_plot('times',16,2,8)
legend off
yyaxis right
plot(X/1e-15,I)
xlabel('time (fs)')
ylabel('Current (A)')
enhance_plot('times',16,2,8)
legend off

t_core_min = 5e-15;
t_core_max = 14e-15;

x_core = x_ecut2((t_ecut2>t_core_min)&(t_ecut2<t_core_max));
xp_core = xp_ecut2((t_ecut2>t_core_min)&(t_ecut2<t_core_max));
y_core = y_ecut2((t_ecut2>t_core_min)&(t_ecut2<t_core_max));
yp_core = yp_ecut2((t_ecut2>t_core_min)&(t_ecut2<t_core_max));
p_core = p_ecut2((t_ecut2>t_core_min)&(t_ecut2<t_core_max));

x_ecut2 = x_ecut2 - mean(x_core);
xp_ecut2 = xp_ecut2 - mean(xp_core);
y_ecut2 = y_ecut2 - mean(y_core);
yp_ecut2 = yp_ecut2 - mean(yp_core);

x_core = x_core - mean(x_core);
xp_core = xp_core - mean(xp_core);
y_core = y_core - mean(y_core);
yp_core = yp_core - mean(yp_core);

gamm0 = mean(p_core);

emitx0=sqrt(mean(x_core.^2).*mean(xp_core.^2)-mean(x_core.*xp_core).^2).*gamma0;
emity0=sqrt(mean(y_core.^2).*mean(yp_core.^2)-mean(y_core.*yp_core).^2).*gamma0;
betax0=mean(x_core.*x_core).*gamma0./emitx0;
betay0=mean(y_core.*y_core).*gamma0./emity0;
alphax0=-mean(x_core.*xp_core).*gamma0./emitx0;
alphay0=-mean(y_core.*yp_core).*gamma0./emity0;


%remove old correlation
xp_ecut2=xp_ecut2+alphax0.*x_ecut2./betax0;
yp_ecut2=yp_ecut2+alphay0.*y_ecut2./betay0;

%scale the beam
x_ecut2=x_ecut2.*sqrt(mbetax./betax0);
y_ecut2=y_ecut2.*sqrt(mbetay./betay0);
xp_ecut2=xp_ecut2.*sqrt(betax0./mbetax);
yp_ecut2=yp_ecut2.*sqrt(betay0./mbetay);

%add new correlation
xp_ecut2=xp_ecut2-malphax.*x_ecut2./mbetax;
yp_ecut2=yp_ecut2-malphay.*y_ecut2./mbetay;

x_core = x_ecut2((t_ecut2>t_core_min)&(t_ecut2<t_core_max));
xp_core = xp_ecut2((t_ecut2>t_core_min)&(t_ecut2<t_core_max));
y_core = y_ecut2((t_ecut2>t_core_min)&(t_ecut2<t_core_max));
yp_core = yp_ecut2((t_ecut2>t_core_min)&(t_ecut2<t_core_max));
p_core = p_ecut2((t_ecut2>t_core_min)&(t_ecut2<t_core_max));


emitx=sqrt(mean(x_core.^2).*mean(xp_core.^2)-mean(x_core.*xp_core).^2).*gamma0;
emity=sqrt(mean(y_core.^2).*mean(yp_core.^2)-mean(y_core.*yp_core).^2).*gamma0;
betax=mean(x_core.*x_core).*gamma0./emitx;
betay=mean(y_core.*y_core).*gamma0./emity;
alphax=-mean(x_core.*xp_core).*gamma0./emitx;
alphay=-mean(y_core.*yp_core).*gamma0./emity;

% DG = 10;
% t_spoil_min = 5e-15;
% t_spoil_max = 12.5e-15;
% spoil_p = DG*randn(size(p_ecut2));
% spoil_p((t_ecut2>t_spoil_min)&(t_ecut2<t_spoil_max)) = 0;
% p_ecut2 = p_ecut2 + spoil_p;



pdata_new = [x_ecut2,xp_ecut2,y_ecut2,yp_ecut2,z,p_ecut2];
nslice = (max(pdata_new(:,5))-min(pdata_new(:,5)))/xlamds;

[X,Y,Z,dx,dy] = contour_plot(t_ecut2,p_ecut2,300,300,0);
I = sum(Z,1)/dx*Q/Np;

figure
yyaxis left
imagesc(X/1e-15,Y,Z)
axis xy
xlabel('time (fs)')
ylabel('\gamma')
colormap(jetvar)
enhance_plot('times',16,2,8)
legend off
yyaxis right
plot(X/1e-15,I)
xlabel('time (fs)')
ylabel('Current (A)')
enhance_plot('times',16,2,8)
legend off


slicepara = slice_beam(x_ecut2,xp_ecut2,y_ecut2,yp_ecut2,t_ecut2,p_ecut2,300,Q/length(t)*length(t_ecut2));

figure
plot(slicepara.tbin/1e-15,slicepara.current)
xlabel('t (fs)')
ylabel('current (A)')
enhance_plot('times',16,2,8);
legend off

figure
yyaxis left
plot(slicepara.tbin/1e-15,slicepara.avg_gamma)
xlabel('t (fs)')
ylabel('average \gamma')
enhance_plot('times',16,2,8);
legend off
yyaxis right
plot(slicepara.tbin/1e-15,slicepara.std_gamma)
xlabel('t (fs)')
ylabel('std \gamma')
enhance_plot('times',16,2,8);
legend off

figure
yyaxis left
plot(slicepara.tbin/1e-15,slicepara.emitx)
xlabel('t (fs)')
ylabel('emittance x')
enhance_plot('times',16,2,8);
legend off
yyaxis right
plot(slicepara.tbin/1e-15,slicepara.emity)
xlabel('t (fs)')
ylabel('emittance y')
enhance_plot('times',16,2,8);
legend off

figure
yyaxis left
plot(slicepara.tbin/1e-15,slicepara.betax)
xlabel('t (fs)')
ylabel('\beta_x')
enhance_plot('times',16,2,8);
legend off
yyaxis right
plot(slicepara.tbin/1e-15,slicepara.betay)
xlabel('t (fs)')
ylabel('\beta_y')
enhance_plot('times',16,2,8);
legend off

figure
yyaxis left
plot(slicepara.tbin/1e-15,slicepara.alphax)
xlabel('t (fs)')
ylabel('\alpha_x')
enhance_plot('times',16,2,8);
legend off
yyaxis right
plot(slicepara.tbin/1e-15,slicepara.alphay)
xlabel('t (fs)')
ylabel('\alpha_y')
enhance_plot('times',16,2,8);
legend off

figure
yyaxis left
plot(slicepara.tbin/1e-15,slicepara.sigx)
xlabel('t (fs)')
ylabel('\sigma_x')
enhance_plot('times',16,2,8);
legend off
yyaxis right
plot(slicepara.tbin/1e-15,slicepara.sigy)
xlabel('t (fs)')
ylabel('\sigma_y')
enhance_plot('times',16,2,8);
legend off

figure
yyaxis left
plot(slicepara.tbin/1e-15,slicepara.xc)
xlabel('t (fs)')
ylabel('xc')
enhance_plot('times',16,2,8);
legend off
yyaxis right
plot(slicepara.tbin/1e-15,slicepara.yc)
xlabel('t (fs)')
ylabel('yc')
enhance_plot('times',16,2,8);
legend off

if dumpQ
    
    foldname = 'genesis/s4';
    if ~exist(foldname,'dir')
        mkdir(foldname);
    end
    
    distfileout = [foldname,'/beam_core_match.dist'];
    
    fid3 = fopen(distfileout,'w');
    fprintf(fid3, '%12s\n', '# Double-Taper Attosection XFEL');
    fprintf(fid3, '%12s\n', '? version = 1.0');
    fprintf(fid3, '%11s  %e\n', '? charge = ',Q/np*length(z));
    fprintf(fid3, '%9s  %d\n', '? size = ',length(z));
    fprintf(fid3, '%12s\n', '? COLUMNS X XPRIME Y YPRIME Z GAMMA');
    
    fprintf(fid3, '%12.10E\t %12.10E\t %12.10E\t %12.10E\t %12.10E\t %12.10E\t\n',pdata_new');
    fclose(fid3);
    
    fclose all;
end


