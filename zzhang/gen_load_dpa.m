function part6d = gen_load_dpa(mainoutfn,dpafn)

if nargin<2||isempty(dpafn)
  dpafn = [mainoutfn,'.dpa'];
end

fid=fopen(dpafn);
npart = gen_get_par(mainoutfn,'npart');
zsep = gen_get_par(mainoutfn,'zsep');
nslice = gen_get_par(mainoutfn,'nslice');
xlamds = gen_get_par(mainoutfn,'xlamds');
Np = npart*nslice;
raw_data = fread(fid,Np*6,'double');
fclose(fid);

dpa = reshape(raw_data,npart,6,nslice);

part6d = zeros(Np,6);
part6d(:,6) = reshape(dpa(:,1,:),Np,1);
gamavg = mean(part6d(:,6));
part6d(:,1) = reshape(dpa(:,3,:),Np,1);
part6d(:,2) = reshape(dpa(:,5,:),Np,1)/gamavg;
part6d(:,3) = reshape(dpa(:,4,:),Np,1);
part6d(:,4) = reshape(dpa(:,6,:),Np,1)/gamavg;

for is = 1:nslice
    part6d(((is-1)*npart+1):(is*npart),5) = (dpa(:,2,is)+2*pi*zsep*(is-1))/(2*pi)*xlamds;
end
% part6d(:,5) = part6d(:,5) - mean(part6d(:,5));