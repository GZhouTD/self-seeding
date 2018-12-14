clc;clear
% close all

ssnum = 28;
foldnum = 1;
plot_energy = 1;
zshow = [0:0.3:37figure];
snum = 51;
R56 = 150;
filter = 3;
distmat = ['genesis/XLBEG2_ssnum_',num2str(snum),'_chicane_',num2str(R56),'_filter_',num2str(filter),'/beam_core_match.mat'];
if ~exist(distmat,'file')
    [pdata,Q,npart] = read_distfile(['genesis/XLBEG2_ssnum_',num2str(snum),'_chicane_',num2str(R56),'_filter_',num2str(filter),'/beam_core_match.dist']);
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
    [data,info] = readOutput(['genesis/XLBEG2_ssnum_',num2str(snum),'_chicane_',num2str(R56),'_filter_',num2str(filter),'/s',num2str(ssnum(ss)),'/genesis.out']);
    paras = load(['genesis/XLBEG2_ssnum_',num2str(snum),'_chicane_',num2str(R56),'_filter_',num2str(filter),'/s',num2str(ssnum(ss)),'/locallog']);
    for ii = 1:length(zshow)
        z = data.z;
        dz = abs(z - zshow(ii));
        power = data.power(dz == min(dz),:)/1e9;
        powermax = max(power);
        t = data.t/3e8/1e-15;
        tmax = max(t);
        ff = data.freq;
        freq0=299792458/info.lambda*6.62606957e-34/1.60217657e-19;
        df = (ff-freq0)/freq0;
        spectrum = data.spectrum(dz == min(dz),:);
        figure(10)
        subplot(2,1,1)
        plot(t,power)
        xlabel('Time (fs)')
        ylabel('Power (GW)')
        set(gca,'xlim',[0, tmax]);
        title(['z = ',num2str(zshow(ii)),'m']);
        enhance_plot('times',16,2,8);
        legend off
        
        subplot(2,1,2)
        farfield_spec = getSpectrum(data.farfield(dz == min(dz),:),data.signalphase(dz == min(dz),:));
        plot(df,farfield_spec)
        set(gca,'xlim',[-0.1,0.1])
        xlabel('\Delta\omega/\omega')
        ylabel('Spectra (a.u.)')
        enhance_plot('times',16,2,8)
        legend off
        
%         print('-dpng','-r72',['figures/z_',num2str(zshow(ii)/0.3+1),'.png']);
        %     spec_center(ii) = sum(farfield_spec.*df)./sum(farfield_spec);
        
        if ii == 0
            pause;
        else
            pause(0.1)
        end
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
set(gca,'xlim',[0,37])
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