%% -------------------------------------------------------------- %%
%                 Spatiotemporal Response function of crystal             %
%% -------------------------------------------------------------- %%
function [GO,GH,t] = Response(d,Ep,theta,chi0)
% ----------------------------------------------------------------------- %
% 
% INPUT:
%  d,     The thickness of crystal, [m]
%  Ep,    photon energy, [eV]
%  theta, Bragg ange, [deg]
%  extin, Extinction length, [m]
%  chi0,  Electric susceptibility of the material
% ----------------------------------------------------------------------- %
% d=30e-6;
% Ep=4000;
% theta=48.924;
% extin=2.8688267611461509e-6;
% chi0 = 0.33863e-04 +1i*0.53162e-06;  %diamond electric susceptibility (111)
% chi0 = 0.66596e-5+1i*0.19046e-6;    %silicon electric susceptibility (111)

h = 6.62607004e-34;
a = 1.602176565e-19;
b = -1;
c = 299792458; 
l = h*c/Ep/a;            % wavelength
w0 = 2*pi*c/l;
K0 = w0/c;
t = linspace(0,300,20000)*1e-15;
extin = sind(theta)/(K0*abs(chi0))
gamma0 = cosd(theta-90);
gammaH = cosd(theta+90);
C = exp(1i*chi0*K0*d/(2*gamma0));
Ta= 2*extin/c*sqrt(abs(b))*sind(theta);
T0 = 2*extin^2/(c*d/gamma0);
Td = 2*d*(sind(theta))^2/(c*abs(gammaH));
wH = -chi0/(2*(sind(theta))^2)*(b-1)/(2*b);
x=sqrt(t/T0.*(1+t/Td));
g = sqrt(abs(b))*abs(chi0)/chi0;
%G = c*K0*chi0/(2*(sin(theta))^2);
gh=1j*g/Ta*(besselj(1,t/Ta)./(t/Ta)).*exp(-1i*wH*w0*t);
GH=gh./max(abs(gh));
go = -C/(2*T0).*(besselj(1,x)./x).*exp(-1i*wH*w0*t);
GO = go./max(abs(gh));







