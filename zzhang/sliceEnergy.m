function [zbin,Eavg,Erms] = sliceEnergy(z,E,Nz)
naprt = length(z);
zmin = min(z);
zmax = max(z);
dz = (zmax - zmin)/Nz;
zhist = (zmin-dz/2):dz:(zmax+dz/2);
Eavg = zeros(length(zhist)-1,1);
Erms = zeros(length(zhist)-1,1);
for ii = 2:length(zhist)
    Ebin = E(z > zhist(ii-1) & z<zhist(ii));
    Eavg(ii-1) = mean(Ebin);
    Erms(ii-1) = std(Ebin);
end

zbin = (zhist(1:end-1)+zhist(2:end))/2;