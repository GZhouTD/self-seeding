clc;clear
% close all

[data,info] = readOutput('genesis/XLBEG2_ssnum_51_chicane_100_filter_3/s2/genesis.out');


%%
% close all
figure
plot(data.z,mean(data.xsize,2),'b')
hold on
plot(data.z,mean(data.ysize,2),'r')
hold off

figure
plot(data.z,data.aw)


figure
semilogy(data.z,mean(data.power,2))


z_monitor = 30;
dz = abs(data.z-z_monitor);
figure
plot(data.freq,data.spectrum(dz == min(dz),:))

freqmax = data.freq(data.spectrum(dz==min(dz),:) == max(data.spectrum(dz==min(dz),:)));

freq = data.freq - freqmax;
spectrum = data.spectrum(dz==min(dz),:);
spec = spectrum./max(spectrum);

power = data.power(dz == min(dz),:);
figure
plot(data.t/3e8/1e-15,power/1e9,'r')
xlabel('t (fs)')
ylabel('Power (GW)')

freq0=299792458/info.lambda*6.62606957e-34/1.60217657e-19;
ff = data.freq;
df = (ff-freq0)/freq0;

farfield_spec = getSpectrum(data.farfield(dz == min(dz),:),data.signalphase(dz == min(dz),:));
figure
plot(df,farfield_spec)

figure
plot(data.t/3e8/1e-15,data.power(dz==min(dz),:))


