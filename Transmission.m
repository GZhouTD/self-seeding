%% -------------------------------------------------------------- %%
%
%% -------------------------------------------------------------- %%
function [R001,R00,R0H] = Transmission(crystal_struc,pho_energy,cry_bragg,freq_arry)

% Ep = 9830;
% dE = 0.4;
% d = 50e-6;
% theta = 45;
format long;
global h_Plank e_charge c_speed

cry_thickness = crystal_struc.cry_thickness;        % crystal thickness [m]
cry_asymmetry = crystal_struc.cry_asymmetry;        % asymmetric angle[deg]
ele_suscept0  = crystal_struc.ele_suscept0;         % electric susceptibility
ele_susceptH  = crystal_struc.ele_susceptH;         % electric susceptibility
ele_susceptHbar  = crystal_struc.ele_susceptHbar;         % electric susceptibility
% chi0 = 0.33863e-04 +1i*0.53162e-06;   %C(111)
% chi0 = 0.38022e-05+1i*0.21732e-07;   %C(004)
% Ep = pho_energy;
% dE = 1000;
% Ome = ang_freq*dE/Ep;
% %DE = linspace(-dE,dE,500000);
% Omega = linspace(-Ome,Ome,500000);

gamma0 = cosd(cry_bragg+cry_asymmetry-90);          % direction cosine
gammaH = cosd(cry_bragg-cry_asymmetry+90);          % direction cosine
asy_factor = gamma0/gammaH;                         % asymmetric factor
wavelength = h_Plank*c_speed/pho_energy/e_charge;   % wavelength [m]
ang_freq = 2*pi*c_speed/wavelength;                 % central angular frequency [Hz]
wave_num = ang_freq/c_speed;                        % wave number


extin_len = sqrt(gamma0*abs(gammaH))/(wave_num*sqrt(ele_susceptH*ele_susceptHbar));  % extinction length
A = cry_thickness/extin_len;
C = exp(1i*ele_suscept0*wave_num*cry_thickness/(2*gamma0));
G = sqrt(abs(asy_factor)*(ele_susceptH*ele_susceptHbar))/ele_susceptHbar;
Omega = 2*pi*freq_arry-ang_freq;                    % frequency deviation
y = wave_num*extin_len/(2*gamma0)*(asy_factor*(-4*Omega*(sind(cry_bragg))^2/ang_freq).*(1-2*Omega/ang_freq)...
    +ele_suscept0*(1-asy_factor));
Y1 = -y-sqrt(y.^2+asy_factor/abs(asy_factor));
Y2 = -y+sqrt(y.^2+asy_factor/abs(asy_factor));
R1 = G*Y1;
R2 = G*Y2;
R00 = exp(1i*(ele_suscept0*wave_num*cry_thickness/2/gamma0+A/2*Y1)).*(R2-R1)./(R2-R1.*exp(1i*A/2*(Y1-Y2)));    % forward diffraction amplitude
R0H = R1.*R2.*(1-exp(1i*A/2*(Y1-Y2)))./(R2-R1.*exp(1i*A/2*(Y1-Y2)));       %   bragg diffraction amplitude
R001 = R00-C;                               % actual forward diffraction
%figure
%plot(DE,abs(R001).^2)
%plotyy(DE,abs(R00),DE,unwrap(angle(R00)))




