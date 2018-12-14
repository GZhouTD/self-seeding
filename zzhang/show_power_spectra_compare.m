clc;clear
close all

ssnum1 = 51:100;
ssnum2 = ssnum1;
zshow = 27.45;
snum = 51;
R56 = 100;
filter = 3;
foldname =  ['genesis/dump_1_XLBEG_ssnum_',num2str(snum),'_chicane_',num2str(R56),'_filter_',num2str(filter)];
for ss = 1:length(ssnum1)
    [data1,info1] = readOutput([foldname,'/s',num2str(ssnum1(ss)),'/genesis.out']);
    [data2,info2] = readOutput([foldname,'/s',num2str(ssnum2(ss)),'/genesis3_15.out']);
    for ii = 1:length(zshow)
        z = data1.z;
        z2 = data2.z;
        dz = abs(z - zshow(ii));
        dz2 = abs(z2 - zshow(ii));
        power1 = data1.power(dz == min(dz),:)/1e9;
        power2 = data2.power(dz2 == min(dz2),:)/1e9;
        powermax1 = max(power1);
        powermax2 = max(power2);
        t1 = data1.t/3e8/1e-15;
        tmax = max(t1);
        t2 = data2.t/3e8/1e-15;
        ff1 = data1.freq;
        ff2 = data2.freq;
        freq1=299792458/info1.lambda*6.62606957e-34/1.60217657e-19;
        freq2=299792458/info2.lambda*6.62606957e-34/1.60217657e-19;
        df1 = (ff1-freq1)/freq1;
        df2 = (ff2-freq2)/freq2;
        
        figure(200)
%         set(gcf,'position',[100,100,800,800])
%         subplot(3,1,1)
        plot(t1,power1)
        hold on
        plot(t2,power2)
        hold off
        xlabel('Time (fs)')
        ylabel('Power (GW)')
        legend('1','2','location','best')
        set(gca,'xlim',[0, tmax]);
        enhance_plot('times',16,2,8);
        %     legend off
        title(['s = ',num2str(ssnum1(ss)),' and ',num2str(ssnum2(ss)),' z = ',num2str(zshow(ii)),'m']);
%         subplot(3,1,2)
%         plot(t1,data1.energy(dz == min(dz),:))
%         hold on
%         plot(t2,data2.energy(dz2 == min(dz2),:))
%         hold off
%         ylabel('energy')
%         enhance_plot('times',16,2,8);
% %         ylim([-100,60])
%         legend off
        
        
%         subplot(3,1,3)
% %         farfield_spec1 = getSpectrum(data1.farfield(dz == min(dz),:),data1.signalphase(dz == min(dz),:));
% %         farfield_spec2 = getSpectrum(data2.farfield(dz == min(dz),:),data2.signalphase(dz2 == min(dz2),:));
%         farfield_spec1 = data1.spectrum(dz == min(dz),:);
%         farfield_spec2 = data2.spectrum(dz2 == min(dz2),:);
%         plot(df1,farfield_spec1)
%         hold on
%         plot(df2,farfield_spec2)
%         hold off
%         set(gca,'xlim',[-0.2,0.2])
%         xlabel('\Delta\omega/\omega')
%         ylabel('auto unit')
%         enhance_plot('times',16,2,8)
%         legend off
        
        if ii == 0
            pause;
        else
            pause(0.1)
        end
    end
end
