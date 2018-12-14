function [dE,zc,sigz] = long_wake(z,L,Q,Nbin,Wz)

%        [dE,zc,sigz] = long_wake(z[,L,Ne,Nbin,fn,pcwake])
%
%	Function to return the wakefield induced energy profile vs. z for
%	a set of given axial coordinates "z".
%
%  INPUTS:	z:		The internal axial coordinates, within the bunch, of
%					each electron with respect to any fixed point [m]
%			L:		(Optional, DEF=1 m) The length of the linac [m]
%			Ne:		(Optional, DEF=1  ) The number of electrons in the bunch
%			Nbin:   (Optional, DEF=100) The number of bins to use
%			fn:		(Optional, DEF='slac.dat') File name containing longitudinal
%					point wake (DEF='slac.dat')
%			pcwake:	(Optional, DEF=none) Point-charge wake used instead of file
%
%  OUTPUTS:	dE:		The energy loss per binned bunch slice [MeV]
%			zc:		The sample points along z where dE is calculated [m]
%					[e.g. plot(zc,dE)]
%			sigz:	rms bunch length (standard deviation) [m]
%
%           (2013-02-25) This version revised by Tim Maxwell primarily
%           per Spencer Gessner's changes to optimize redundant nested for-loop.
%           Modifications to preload wakefield externally not included.
%           Forced extrapolation of wake to zero added, with warning.
%=============================================================================

% Number of simulated particles.
nn   = length(z);
% RMS z-spread
sigz = std(z);

if nn < 100
  disp(' ')
  disp('Probably too few particles in your input "z" array')            
  disp(' ')
end
if nn > 15E6
  disp(' ')
  disp('>5E6 particles in your "z" array - let''s not get carried away now')
  disp(' ')
end
if sigz < 1E-6
  disp(' ')
  disp('Bunch length of <1 micron may be too short for this Green''s function')
  disp(' ')
end

if ~exist('L')
  L = 1; 					% default length of S-band linac [m]
end
if ~exist('Ne')
  Q = 250e-12; 					% default number of e- in bunch      
end
if ~exist('Nbin')
  Nbin = 300;  				% default number simulation particles
end

zfnvar = Wz(:,1)*3e8;			% m
Wfnvar = Wz(:,2);			% V/C/m
% zc = zfnvar;
nA     = length(zfnvar);
% Histogram particles into Nbins with zeros at ends for interpolation
zmin = min(z);
zmax = max(z);
dzc = (zmax - zmin)/Nbin;
zc = (zmin-dzc/2):dzc:(zmax+dzc/2);
Nbin = length(zc);
N = hist(z,zc);
dzc = zc(2)-zc(1);			% (TJM 2013-02-25) Above line replaced: hist returns constant-spacing in zc.  (Credit: Spencer)

% Add zero padding to close ends. (Note that extrapolation beyond limits is
% now explicitly set to zero below. using EXTRAPVAL parameter for interp1.)

% maxz_fn = zfnvar(nA-1);			% max Z in wake file (last point usually projected estimate)
% % if (max(z)-min(z)) > maxz_fn
% if (zc(Nbin)-zc(1)) > maxz_fn % (TJM 2013-02-25) Above line replaced: Bunch extents already defined by binned axis. (Credit: Spencer)
%   disp(' ')
%   if ~exist('pcwake')
%     disp(['WARNING: maximum axial spread is > ' num2str(maxz_fn*1e3) ' mm and ' fn ' is inaccurate there.'])
%     disp('         Wake will be extrapolated to zero as needed');
%   else
%     disp(['WARNING: maximum axial spread is > ' num2str(maxz_fn*1e3) ' mm and the RW-wake is inaccurate there'])
%     disp('         Wake will be extrapolated to zero as needed');
%   end
%   disp(' ')
% end

% delta vector
dE  = zeros(Nbin,1);
% electron charge (SI units, coulombs)
e   = 1.602176565E-19;
% Beam current normalization factor times drift length.

scl = (Q/nn)*L;

% Following code replaces commented block above.
N = N.';
zc = zc.';
% Bin separation vector
dzi = zc - max(zc);
% Interpolated wake vector
Wf = interp1(zfnvar,Wfnvar,dzi,'linear',0);
% Self-wake bin
Wf(1) = Wf(1)/2;
% Sum delta due to wake from leading bins
figure
plot(zc/3e8/1e-15,N*scl)
for j =Nbin:-1:1
    dE(j) = sum(scl*N(j:1:Nbin).*Wf((Nbin-j+1):-1:1));
end
% End of new block.
