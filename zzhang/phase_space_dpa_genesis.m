clc;clear
close all

filename = 'genesis/XLBEG2_ssnum_51_chicane_100_filter_3/s1/genesis.out';
dpa = gen_load_dpa(filename);
[data,info] = readOutput(filename);

R56 = (0:10:200)*1e-6;
T566 = -1.5*R56;
clight = 2.998e8;
t = dpa(:,5);
p = dpa(:,6);
delta = (p - mean(p))/mean(p);
nbin = 300;
Is = zeros(length(R56),nbin);
Imax = zeros(size(R56));
tbin = zeros(length(R56),nbin);

for ii = 1:length(R56)
    ii
    t2 = t + delta*R56(ii);
    [I,tt] = current_plot_dpa(t2/clight,data,info,nbin);
%     figure
%     plot(tt/1e-15,I)
%     xlabel('t (fs)')
%     ylabel('current (A)')
%     title(['R56 = ',num2str(R56(ii)*1e6),' \mum']);
%     enhance_plot('times',16,2,8);
%     legend off
    Is(ii,:) = I;
    tbin(ii,:) = tt;
    Imax(ii) = max(I(abs(tt)<2e-15));
end

figure
plot(R56/1e-6,Imax/1000)
xlabel('R56 (\mum)')
ylabel('Current peak (kA)')
enhance_plot('times',16,2,8)
legend off
