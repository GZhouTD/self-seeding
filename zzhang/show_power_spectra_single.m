clc;clear
close all
ssnum = 9;
% ssnum = 51:60;
% ssnum = 171:175;
% ssnum = 91:100;
% ssnum = 71:80;
foldnum = 1;
plot_energy = 1;
nu_stop = 8;
% zshow = [7.35:0.15:27.75];
zshow = nu_stop*3.3+0.3+0.15*(nu_stop-1);
% zshow = 38.1;
% zshow = 45;
% zshow = 24.3;
gamma_modulation = 50;
taper = 10e-4;
snum = 51;
R56 = 100;
filter = 3;
foldname =  ['genesis/dump_1_XLBEG_ssnum_',num2str(snum),'_chicane_',num2str(R56),'_filter_',num2str(filter)];

distmat = [foldname,'/beam_core_match.mat'];
if ~exist(distmat,'file')
    [pdata,Q,npart] = read_distfile([foldname,'/beam_core_match.dist']);
    z = pdata{5};
    z = z - min(z);
    gamma = pdata{6};
    [X,Y,Z,dx,dy] = contour_plot(z/3e8*1e15,gamma,300,300,0,Q);
    I = sum(Z,1)*Q/npart/(dx*1e-15);
    save(distmat,'X','Y','Z','dx','dy','I','Q','npart');
else
    load(distmat);
end

for ss = 1:length(ssnum)
[data,info] = readOutput([foldname,'/s',num2str(ssnum(ss)),'/genesis3_15.out']);
% [data,info] = readOutput([foldname,'/s',num2str(ssnum(ss)),'/genesis2_',num2str(taper),'.out']);
% [data,info] = readOutputharmBPFFP([foldname,'/s',num2str(ssnum(ss)),'/genesis.out']);
paras = load([foldname,'/s',num2str(ssnum(ss)),'/locallog']);
% paras(12)
for ii = 1:length(zshow)
    z = data.z;
    dz = abs(z - zshow(ii)); 
    power = data.power(dz == min(dz),:)/1e9;
    powermax = max(power);
%     [pks,locs] = findpeaks(power,'MinPeakHeight',powermax/5,'NPeaks',5);
    t = data.t/3e8/1e-15;
%     t_peaks = t(locs);
    tmax = max(t);
    ff = data.freq;
    freq0=299792458/info.lambda*6.62606957e-34/1.60217657e-19;
    df = (ff-freq0)/freq0;
    spectrum = data.spectrum(dz == min(dz),:);
    figure(100)
    set(gcf,'position',[100,100,800,800])
    subplot(3,1,1)
    imagesc(X,Y,Z);
    axis xy
    xlabel('Time (fs)')
    ylabel('\gamma');
    set(gca,'xlim',[0, tmax]);
%     title(['s = ',num2str(ssnum(ss)),' z = ',num2str(zshow(ii)),'m taper1 = ',num2str(paras(10)),' taper2 = ',num2str(paras(25))]);
    title(['s = ',num2str(ssnum(ss)),' z = ',num2str(zshow(ii)),'m taper1 = ',num2str(paras(5))]);
    enhance_plot('times',16,2,8);
    legend off
    colormap(jetvar);
    subplot(3,1,2)
    yyaxis left
    plot(t,power)
    xlabel('Time (fs)')
    ylabel('Power (GW)')
    set(gca,'xlim',[0, tmax]);
%     title(['s = ',num2str(ssnum(ss)),' z = ',num2str(zshow(ii)),'m t12 = ',num2str(paras(5)),' t5 = ',num2str(paras(9))]);
    
    title(['s = ',num2str(ssnum(ss)),' z = ',num2str(zshow(ii)),'m taper1 = ',num2str(paras(5))]);
    enhance_plot('times',16,2,8);
    legend off
    if plot_energy
        yyaxis right
        plot(t,data.energy(1,:),'b--')
        hold on
        plot(t,data.energy(dz == min(dz),:))
        hold off
        ylabel('energy')
        enhance_plot('times',16,2,8);
        %     ylim([-100,60])
        legend off
    else
        yyaxis right
        plot(t,data.bunching(dz == min(dz),:))
        ylabel('bunching')
        enhance_plot('times',16,2,8);
        %     ylim([-100,60])
        legend off
    end

    subplot(3,1,3)
    farfield_spec = getSpectrum(data.signal(dz == min(dz),:),data.signalphase(dz == min(dz),:));
    ffspectrum = farfield_spec;
%     ffspectrum = data.ffspectrum(dz == min(dz),:);
    plot(df,ffspectrum)
    set(gca,'xlim',[-0.2,0.2])
    xlabel('\Delta\omega/\omega')
    ylabel('auto unit')
    enhance_plot('times',16,2,8)
    legend off
    
%     print('-dpng','-r72',['ideal/s',num2str(foldnum),'/figures/result_',num2str(ssnum(ss)),'.png']);
    spec_center(ii) = sum(farfield_spec.*df)./sum(farfield_spec);
    
    if ii == 0
        pause;
    else
        pause(0.1)
    end
end
% figure
% plot(data.z,mean(data.xsize,2))
% hold on
% plot(data.z,mean(data.ysize,2))
% xlabel('z (m)')
% ylabel('beam size')
% enhance_plot('times',16,2,8)
% legend off
% 
% figure
% plot(data.z,mean(data.power,2))
% xlabel('z (m')
% ylabel('Power (W)')
% enhance_plot('times',16,2,8)
% legend off
% 
% figure
% plot(data.t/3e8/1e-15,data.energy)
% xlabel('t (fs)')
% ylabel('energy')
% enhance_plot('times',16,2,8)
% legend off
% 
% 
aw = data.aw;
zz = data.z;
zz(aw<1) = [];
aw(aw<1) = [];
daw = diff(aw);
dzz = diff(zz);
z_avg = (zz(2:end) + zz(1:(end-1)))/2;
daw_dz = daw./dzz/2.5;
figure(101)
% yyaxis left
plot(zz,aw,'.')
xlabel('z (m)')
ylabel('aw')
% set(gca,'ylim',[2.34,2.4])
enhance_plot('times',16,2,8)
legend off
% yyaxis right
% plot(z_avg,daw_dz)
% ylabel('taper')
% enhance_plot('times',16,2,8)
% legend off

% figure(102)
% plot(zshow,spec_center)
% xlabel('z (m)')
% ylabel('spectra center')
% % set(gca,'ylim',[-0.02,0.02])
% enhance_plot('times',16,2,8)
% legend off

end