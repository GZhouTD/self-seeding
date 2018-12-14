clc;clear
close all

ssnum = 51;
R56 = 100;
filter = 3;

foldname = ['elegant/scan/BC1COL/dump_10M/col_3_s',num2str(ssnum)];

filename = [foldname,'/XLBEG_chicane_',num2str(R56),'_filter_',num2str(filter),'_match2_.wig'];

pdata = load(filename);
Q = 1.8770e-10;

%%
x = pdata(:,1);
xp = pdata(:,2);
y = pdata(:,3);
yp = pdata(:,4);
t = pdata(:,5);
p = pdata(:,6);
clear pdata

t = t - mean(t);

t = -t;

[X,Y,Z,~,~,I] = contour_plot_current(t,p*0.511,300,300,0,Q);
figure
set(gcf,'position',[100,100,800,800])
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
enhance_plot('times',16,2,8)
legend off


pcut1 = 5010/0.511;
tcut1 = 0;
tcut2 = 20e-15;

tbake = t;
pbake = p;

x_cut = x((tbake>tcut1)&(tbake<tcut2)&(pbake<pcut1));
xp_cut = xp((tbake>tcut1)&(tbake<tcut2)&(pbake<pcut1));
y_cut = y((tbake>tcut1)&(tbake<tcut2)&(pbake<pcut1));
yp_cut = yp((tbake>tcut1)&(tbake<tcut2)&(pbake<pcut1));
t_cut = t((tbake>tcut1)&(tbake<tcut2)&(pbake<pcut1));
p_cut = p((tbake>tcut1)&(tbake<tcut2)&(pbake<pcut1));

% [X,Y,Z,~,~,I] = contour_plot_current(t_cut,p_cut*0.511,300,300,0,Q/length(t)*length(t_cut));
% figure
% set(gcf,'position',[100,100,800,800])
% subplot(2,1,1)
% imagesc(X*1e15,Y/1000,Z)
% axis xy
% colormap(jetvar)
% xlabel('time (fs)')
% ylabel('E (GeV)')
% enhance_plot('times',16,2,8)
% legend off
% subplot(2,1,2)
% plot(X*1e15,I)
% xlabel('time (fs)')
% ylabel('Current (kA)')
% enhance_plot('times',16,2,8)
% legend off

clear tbake pbake

%%
t_cut = t_cut - min(t_cut);
tspoil = 12.5e-15;
dp = zeros(size(p_cut));
sigp = 50;
dp = sigp*randn(size(dp));
dp(t_cut<tspoil) = 0;
p_cut = p_cut + dp;

[X,Y,Z,~,~,I] = contour_plot_current(t_cut,p_cut*0.511,300,300,0,Q/length(t)*length(t_cut));
figure
set(gcf,'position',[100,100,800,800])
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
enhance_plot('times',16,2,8)
legend off






sampratio = 2;
x_cut = x_cut(1:sampratio:end);
xp_cut = xp_cut(1:sampratio:end);
y_cut = y_cut(1:sampratio:end);
yp_cut = yp_cut(1:sampratio:end);
t_cut = t_cut(1:sampratio:end);
p_cut = p_cut(1:sampratio:end);

% [X,Y,Z,~,~,I] = contour_plot_current(t_cut,p_cut*0.511,300,300,0,Q/length(t)*length(t_cut)*sampratio);
% figure
% set(gcf,'position',[100,100,800,800])
% subplot(2,1,1)
% imagesc(X*1e15,Y/1000,Z)
% axis xy
% colormap(jetvar)
% xlabel('time (fs)')
% ylabel('E (GeV)')
% enhance_plot('times',16,2,8)
% legend off
% subplot(2,1,2)
% plot(X*1e15,I)
% xlabel('time (fs)')
% ylabel('Current (kA)')
% enhance_plot('times',16,2,8)
% legend off


[X,Y,Z,~,~,I] = contour_plot_current(t_cut,xp_cut,300,300,0,Q/length(t)*length(t_cut));
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


clight = 299792458;
pdata_new = [x_cut,xp_cut,y_cut,yp_cut,t_cut*clight,p_cut];

foldgenesis = ['genesis/s2_XLBEG_ssnum_',num2str(ssnum),'_chicane_',num2str(R56),'_filter_',num2str(filter)];
if ~exist(foldgenesis,'dir')
    mkdir(foldgenesis);
end
% slicepara = slice_beam(x_cut,xp_cut,y_cut,yp_cut,t_cut,p_cut,300,Q);
distfileout = [foldgenesis,'/beam_core_match.dist'];

fid3 = fopen(distfileout,'w');
fprintf(fid3, '%12s\n', '# Double-Taper Attosection XFEL');
fprintf(fid3, '%12s\n', '? version = 1.0');
fprintf(fid3, '%11s  %e\n', '? charge = ',Q/length(t)*length(t_cut)*sampratio);
fprintf(fid3, '%9s  %d\n', '? size = ',length(t_cut));
fprintf(fid3, '%12s\n', '? COLUMNS X XPRIME Y YPRIME Z GAMMA');

fprintf(fid3, '%12.10E\t %12.10E\t %12.10E\t %12.10E\t %12.10E\t %12.10E\t\n',pdata_new');
fclose(fid3);

matdata = [foldgenesis,'/matdata.mat'];
save(matdata,'tspoil','sigp','pcut1','tcut1','tcut2');

fclose all;