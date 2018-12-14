clc;clear
close all

ssnum = 1:20;
foldnum = 1;
plot_energy = 1;
zshow = 24.3;
gamma_modulation = 50;
snum = 51;
R56 = 100;
filter = 3;
foldname =  ['genesis/XLBEG2_ssnum_',num2str(snum),'_chicane_',num2str(R56),'_filter_',num2str(filter)];
power_all = zeros(length(ssnum),5000);
amplitude = zeros(length(ssnum),5000);
phase = zeros(length(ssnum),5000);
spectrum_all = zeros(length(ssnum),5000);

for ss = 1:length(ssnum)
    ss
    [data,info] = readOutput([foldname,'/s',num2str(ssnum(ss)),'/genesis.out']);
    for ii = 1:length(zshow)
        z = data.z;
        dz = abs(z - zshow(ii));
        power = data.power(dz == min(dz),:)/1e9;
        powermax = max(power);x
        t = data.t/3e8/1e-15;
        tmax = max(t);
        ff = data.freq;
        freq0=299792458/info.lambda*6.62606957e-34/1.60217657e-19;
        df = (ff-freq0)/freq0;
        spectrum_all(ss,:) = data.spectrum(dz == min(dz),:);
        power_all(ss,:) = power;
        amplitude(ss,:) = data.signal(dz == min(dz),:);
        phase(ss,:) = data.signalphase(dz == min(dz),:);
    end
end