clc;clear
close all

load('scan_R56.mat');
for ii = 1:3:length(R56)
    figure
    plot(tbin(ii,:)/1e-15,Is(ii,:)/1000)
    xlabel('t (fs)')
    ylabel('current (kA)')
    set(gca,'xlim',[-15,15]);
    title(['R56 = ',num2str(R56(ii)*1e6),' \mum']);
    enhance_plot('times',16,2,8);
    legend off
end

figure
plot(R56/1e-6,Imax/1000)
xlabel('R56 (\mum)')
ylabel('Current peak (kA)')
enhance_plot('times',16,2,8)
legend off
