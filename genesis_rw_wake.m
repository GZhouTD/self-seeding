function genesis_rw_wake(shiftDC)
%	to calculate rw_wake for LCLS undulator chamber, Al, rectangle, 5mm gap
%	are hardcoded. 
%   a current profile [z(m) current (A)] will be read, with bunch head on the left (from originally elegant2current used for Genwake by Sven).
%   output will be saved in outputFile.
% shiftDC=1 means to remove the offset, which could be tapered on the real
% machine.
% genesis_rw_wake('undcur.dat','LCLSwake.dat')

currentFile='undbeg-current.dat';
outputFile='wake.beam';
sig  = 3.5e7;  % Al: 'Conductivity (ohm-1*m-1)'
tau  = 8e-15;      % Al: relaxation time
rf   =  1;              % rf=1: rectangle chamber: rf=0: round chamber    
r    =2.5;             % mm, chamber radius 

c  = 2.99792458E8;
Z0 = 120*pi;

[zs Ipk] = textread(currentFile,'%f %f','delimiter',' ');
Q=integrate(zs/c,Ipk);
r  = r*1E-3;
s0 = (2*r^2/(Z0*sig))^(1/3);

f = Ipk/integrate(zs,Ipk);

s = zs - zs(1);
w = rw_wakefield(s,r,s0,tau,rf);

n = length(s);
E = zeros(n,n);
for j = 1:n
  for i = 1:n
    if i==j
      break
    else
      E(i,j) = w(j-i)*f(i);
    end
  end
end

dz = mean(diff(zs));
Ez = Q*sum(E)*dz; % eV/m/

Ez_mean = integrate(zs,f'.*Ez)
Ez_rms  = sqrt(integrate(zs,f'.*(Ez-Ez_mean).^2));
%Ez_rmsg = 100*rw_esprd(E0/1E9,Ne/1E10,L,r,sigz*1E6,sig);

zs=max(zs)-zs;
zs=flipud(zs);
Ez=flipud(Ez');
Ipk=flipud(Ipk);

if (shiftDC==1)
    Ez = Ez-Ez_mean;
end


%tstr = ['AC Resistive-Wall Wake ({\it\tau} = ' sprintf('%4.1f',tau*1E15) ' fs, {\it\sigma_c} = ' sprintf('%4.2f',sig/1E7) '\times10^7' ...
%         ' /\Omega/m, {\itr} = ' sprintf('%4.1f',r*1E3) ' mm'];

out=[zs Ipk Ez];
 header1 = '? VERSION=1.0';
 header2 = ['? SIZE=' num2str(length(zs))];
 header3 = '? COLUMNS ZPOS CURPEAK ELOSS';
 
fid = fopen(outputFile,'wt');
fprintf(fid,'%s\n',header1);
fprintf(fid,'%s\n',header2);
fprintf(fid,'%s\n',header3);
fprintf(fid,'%14.8e %14.8e %14.8e\n',out');
fclose(fid)


figure(11)
plot(zs*1E6,Ez)
xlabel('{\itz} (\mum)'); ylabel('Eloss (eV/m)')
figure(12)
plot(zs*1E6,Ipk/1E3)
xlabel('\mum'); ylabel('kA');

if (shiftDC==1)
figure(13)
plot(zs*1E6,Ez+Ez_mean)
xlabel('{\itz} (\mum)'); ylabel('Eloss (eV/m)')
figure(12)
plot(zs*1E6,Ipk/1E3)
xlabel('\mum'); ylabel('kA');
end


function s = integrate(x,y,x1,x2)

%       s = integrate(x,y[,x1,x2]);
%
%       Approximate the integral of the function y(x) over x from x1 to 
%       x2.  The limits of integration are given by the optional inputs
%       x1 and x2.  If they are not given the integration will range from
%       x(1) to x(length(x)) (i.e. the whole range of the vector x).
%
%     INPUTS:   x:      The variable to integrate over (row or column
%                       vector of sequential data points)
%               y:      The function to integrate (row or column vector)
%               x1:     (Optional,DEF=x(1)) The integration starting point
%               x2:     (Optional,DEF=x(n)) The integration ending point

%===============================================================================

if any(diff(x)<0);
  error('x must be sequentially ordered data')
end
  
x = x(:);
y = y(:);

[ny,cy] = size(y);
[nx,cx] = size(x);

if (cx > 1) | (cy > 1)
  error('INTEGRATE only works for vectors')
end
if nx ~= ny
  error('Vectors must be the same length')
end

if ~exist('x2')
  i2 = nx;
  if ~exist('x1')
    i1 = 1;
  else
    [dum,i1] = min(abs(x-x1));
  end
else
  [dum,i1] = min(abs(x-x1));
  [dum,i2] = min(abs(x-x2));
end

dx = diff(x(i1:i2));
s = sum(dx.*y(i1:(i2-1)));


