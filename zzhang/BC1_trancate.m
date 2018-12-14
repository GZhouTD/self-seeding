clc;clear
close all

[X,Y,Z,~,~] = sddsplotxy('elegant/scan/BC1MID/s1','BC1MID',250e-12);

figure;
imagesc(X,Y,Z)
axis xy
xlabel('x (m)')
ylabel('y (m)')
colormap(jetvar)
enhance_plot('times',16,2,8)
legend off