clc;clear
close all
fid = fopen('elegant/s1_1M/pdata');
np = fscanf(fid,'%d',1);

pdata = fscanf(fid,'%f',np*6);
fclose(fid);
ngroup = np/6;
sq = (1:np)*6-6;
x = pdata(sq+1);
xp = pdata(sq+2);
y = pdata(sq+3);
yp = pdata(sq+4);
t = pdata(sq+5);
p = pdata(sq+6);
Np = length(x);

figure;
contour_plot(t/1e-15,x,300,300,1);