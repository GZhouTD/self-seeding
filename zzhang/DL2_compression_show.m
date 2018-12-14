clc;clear
close all


Q = 250e-12;


fid = fopen('elegant/DL2_200/pdata');
np = fscanf(fid,'%d',1);
pdata = fscanf(fid,'%f',np*6);
fclose(fid);
ngroup = np/6;
sq = (1:np)*6-6;
t = pdata(sq+5);
p = pdata(sq+6);
Np = length(t);
t = t - mean(t);
[X,Y,Z,dx,dy,I1] = contour_plot_current(t,p,300,300,1,Q);


fid = fopen('elegant/DL2_220/pdata');
np = fscanf(fid,'%d',1);
pdata = fscanf(fid,'%f',np*6);
fclose(fid);
ngroup = np/6;
sq = (1:np)*6-6;
t = pdata(sq+5);
p = pdata(sq+6);
Np = length(t);
t = t - mean(t);
[X,Y,Z,dx,dy,I2] = contour_plot_current(t,p,300,300,1,Q);


fid = fopen('elegant/DL2_280/pdata');
np = fscanf(fid,'%d',1);
pdata = fscanf(fid,'%f',np*6);
fclose(fid);
ngroup = np/6;
sq = (1:np)*6-6;
t = pdata(sq+5);
p = pdata(sq+6);
Np = length(t);
t = t - mean(t);
[X,Y,Z,dx,dy,I3] = contour_plot_current(t,p,300,300,1,Q);

fid = fopen('elegant/DL2_420/pdata');
np = fscanf(fid,'%d',1);
pdata = fscanf(fid,'%f',np*6);
fclose(fid);
ngroup = np/6;
sq = (1:np)*6-6;
t = pdata(sq+5);
p = pdata(sq+6);
Np = length(t);
t = t - mean(t);
[X,Y,Z,dx,dy,I4] = contour_plot_current(t,p,300,300,1,Q);


figure
plot(X*1e15,I1)
hold on
plot(X*1e15,I2)
plot(X*1e15,I3)
plot(X*1e15,I4)
hold off
xlabel('time (fs)')
ylabel('Current (kA)')
set(gca,'xlim',[-40 40]);
legend('R56 = 200\mum','R56 = 220\mum','R56 = 280\mum','R56 = 360\mum','location','best');
enhance_plot('times',16,2,8)