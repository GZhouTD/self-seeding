function [cmplxfld,result] = gen_load_dfl(mainoutfn,dflfn)
%**[cmplxfld,result] = gen_load_dfl(mainoutfn,dflfn)
%  To load the radiation field from dfl file.
%<<Input
%  mainoutfn 
%    - Genesis main output filename, from which Genesis simulations
%      parameters are obtained
%  dflfn 
%    - Genesis output dfl filename
%>>Output
%  cmplxfld 
%    - the complex field data
%      ncar x ncar x nslice matrix
%  result
%    - field data for plot / analysis - needs to check later !!!
%    .dpow:     power over the grid area, ncar x ncar x nslice martix
%    .pow:      power (integration of dpow over all trans. grids), 1 x nslice array
%    .powdens:  power density, ncar x ncar x nslice martix
%    .xv:       x coordinates (m), 1 x ncar array
%    .yv:       y coordinates (m), 1 x ncar array
%    .sv:       s coordinates (m), 1 x nslice array
%    .tv:       t coordinates (s), 1 x nslice array
%    .maxpowdens: maximum power density, 1 x nslice array
%    .ctrpowdens: on-axis power density, 1 x nslice array 
%
%**Modification log
%  established on 2015-01-22-Thu (S.Huang)

if nargin<2||isempty(dflfn)
  dflfn = [mainoutfn,'.dfl'];
end

c = 2.99792458E8;
nslice = gen_get_par(mainoutfn,'nslice');
ncar   = gen_get_par(mainoutfn,'ncar');
slicesize = 2*ncar^2; 
dgrid  = gen_get_par(mainoutfn,'dgrid');
ddgrid = 2*dgrid/(ncar-1);
zsep   = gen_get_par(mainoutfn,'zsep');
xlamds = gen_get_par(mainoutfn,'xlamds');
ds = xlamds*zsep;
sv = (0:(nslice-1))*ds;
tv = sv/c;

fid=fopen(dflfn);
cmplxfld=zeros(ncar,ncar,nslice);
for ii = 1:nslice
  data = fread(fid,slicesize,'double');
  re = reshape(data(1:2:slicesize),ncar,ncar)';
  im = reshape(data(2:2:slicesize),ncar,ncar)';
  cmplxfld(:,:,ii)=complex(re,im);
end
fclose(fid);

% needs to convert to real physical quantity later
dpow = cmplxfld.*conj(cmplxfld); % W
pow = reshape(sum(sum(dpow)),1,nslice); % W
powdens = dpow/ddgrid^2; % W/m^2

result.dpow = dpow;
result.pow = pow;
result.powdens = powdens;
result.maxpowdens = reshape(max(max(powdens)),1,nslice);
result.ctrpowdens = reshape(powdens((ncar+1)/2,(ncar+1)/2,:),1,nslice);
result.xv = -dgrid:ddgrid:dgrid;
result.yv = -dgrid:ddgrid:dgrid;
result.sv = sv;
result.tv = tv;
end
