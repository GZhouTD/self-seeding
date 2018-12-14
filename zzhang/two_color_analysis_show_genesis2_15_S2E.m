clc;clear
close all

load('two_color_genesis2_15.mat');
xlambda = 1.1698e-09*(1-0.035)*(1+0.05)*0.9; %% negative taper aw0 = 2.41
badshot = [19,21,31,36,39,40,42,52,57,66,76,82,98];
goodshot = [16,24,29,32,34,37,43,50];

shotnum = 1:length(data);
shotnum(badshot) = [];

fwhm1 = zeros(1,length(shotnum));
tc1 = zeros(1,length(shotnum));
tp1 = zeros(1,length(shotnum));
peak_power1 = zeros(1,length(shotnum));

fwhm2 = zeros(1,length(shotnum));
tc2 = zeros(1,length(shotnum));
tp2 = zeros(1,length(shotnum));
peak_power2 = zeros(1,length(shotnum));


for ii = 1:length(shotnum)
    result = data{shotnum(ii)};
    
    [fwhm1(ii),tc1(ii),tp1(ii),peak_power1(ii)] = single_spike_analysis(result.t1,result.power1);
    [sp_fwhm1(ii),spc1(ii),spp1(ii),sp_power1(ii)] = single_spike_analysis(result.ff1,result.spec1);
    [fwhm2(ii),tc2(ii),tp2(ii),peak_power2(ii)] = single_spike_analysis(result.t2,result.power2);
    [sp_fwhm2(ii),spc2(ii),spp2(ii),sp_power2(ii)] = single_spike_analysis(result.ff2,result.spec2);
end

fwhm1_avg = mean(fwhm1);
fwhm1_std = std(fwhm1);
fwhm2_avg = mean(fwhm2);
fwhm2_std = std(fwhm2);

dtc = abs(tc2 - tc1)-(110*8)*xlambda/3e8/1e-15;
dtp = abs(tp2 - tp1);
dtc_avg = mean(dtc);
dtc_std = std(dtc);

power1_avg = mean(peak_power1);
power1_std = std(peak_power1);
power2_avg = mean(peak_power2);
power2_std = std(peak_power2);

dspc = abs(spc2 - spc1);
pc = polyfit(dspc,dtc,1);
dtc_unc = dtc - (pc(1)*dspc+pc(2));
std(dtc_unc)

% figure
% plot(fwhm1)
% hold on
% plot(fwhm2)
% hold off
% xlabel('shot #')
% ylabel('fwhm (fs)')
% legend(['900 eV, fwhm=',num2str(fwhm1_avg,'%1.2f'),'\pm',num2str(fwhm1_std,'%1.2f'),'fs'],...
%     ['1100 eV, fwhm=',num2str(fwhm2_avg,'%1.2f'),'\pm',num2str(fwhm2_std,'%1.2f'),'fs'])
% enhance_plot('times',16,2,8)
% 
% figure
% plot(dtc,'o-')
% xlabel('shot #')
% ylabel('{\it\Deltat}')
% legend(['{\it\Deltat}=',num2str(dtc_avg,'%1.2f'),'\pm',num2str(dtc_std,'%1.2f'),'fs'])
% enhance_plot('times',16,2,8)
% 
% figure
% plot(peak_power1)
% hold on
% plot(peak_power2)
% hold off
% xlabel('shot #')
% ylabel('fwhm (fs)')
% legend(['900 eV, Power=',num2str(power1_avg,'%1.2f'),'\pm',num2str(power1_std,'%1.2f'),'GW'],...
%     ['1100 eV, Power=',num2str(power2_avg,'%1.2f'),'\pm',num2str(power2_std,'%1.2f'),'GW'])
% enhance_plot('times',16,2,8)
% 
% figure
% plot(dspc,dtc,'.')
% xlabel('{\it\Deltah\nu}  (eV)')
% ylabel('{\it\Deltat}')
% enhance_plot('times',16,2,8)
% legend off

fwhm1_bin = 0.2:0.05:0.9;
fwhm2_bin = 0.2:0.025:0.6;
figure
subplot(2,1,1)
histogram(fwhm1,fwhm1_bin)
xlabel('fwhm (fs)')
ylabel('count')
legend(['900 eV, fwhm=',num2str(fwhm1_avg,'%1.2f'),'\pm',num2str(fwhm1_std,'%1.2f'),'fs'])
enhance_plot('times',16,2,8)
subplot(2,1,2)
histogram(fwhm2,fwhm2_bin)
xlabel('fwhm (fs)')
ylabel('count')
legend(['1100 eV, fwhm=',num2str(fwhm2_avg,'%1.2f'),'\pm',num2str(fwhm2_std,'%1.2f'),'fs'])
enhance_plot('times',16,2,8)

power1_bin = 80:10:260;
power2_bin = 30:10:200;
figure
subplot(2,1,1)
histogram(peak_power1,power1_bin)
xlabel('Power (GW)')
ylabel('count')
legend(['900 eV, Power=',num2str(power1_avg,'%1.2f'),'\pm',num2str(power1_std,'%1.2f'),'GW'],'location','best')
enhance_plot('times',16,2,8)
subplot(2,1,2)
histogram(peak_power2,power2_bin)
xlabel('Power (GW)')
ylabel('count')
legend(['1100 eV, Power=',num2str(power2_avg,'%1.2f'),'\pm',num2str(power2_std,'%1.2f'),'GW'],'location','best')
enhance_plot('times',16,2,8)

dtc_bin = 4.2:0.05:5.3;
figure
histogram(dtc,dtc_bin)
set(gca,'ylim',[0,15])
xlabel('{\it\Deltat} (fs)')
ylabel('count')
legend(['{\it\Deltat}=',num2str(dtc_avg,'%1.2f'),'\pm',num2str(dtc_std,'%1.2f'),'fs'],'location','best')
enhance_plot('times',16,2,8)

shotnum = 1:length(data);
for ii = 1:length(goodshot)
    result = data{shotnum(goodshot(ii))};
    t1 = result.t1+(110*8)*xlambda/3e8/1e-15;
    sp1(ii,:) = result.power1(t1<10.1);
    t2 = result.t2;
    sp2(ii,:) = result.power2(t2>9.9);
    
end
t1 = t1(t1<10.1);
t2 = t2(t2>9.9);
sp1 = sp1/200;
sp2 = sp2/200;
figure;
plot(t1,sp1)
hold on
plot(t2,sp2)

figure
plot(t1,sp1(1,:),'b')
hold on
plot(t2,sp2(1,:),'r')
for ii = 2:length(sp2(:,1))
    plot(t1,sp1(ii,:)+ii-1,'b')
    plot(t2,sp2(ii,:)+ii-1,'r')
end
hold off
set(gca,'xlim',[5,14])
xlabel('time (fs)')
ylabel('Power profile')
enhance_plot('times',16,2,8)
legend off


