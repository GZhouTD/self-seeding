clc;clear
close all
filename = 'wiggler_wake_6per_2um';
pdata = load(['wake/',filename,'.txt']);
t = pdata(:,1);
Wz = pdata(:,2);
% t = t-t(1);
t(1) = 0;
%print longitudinal wake
fn=strcat(['wake/',filename,'.sdds']);
fid = fopen(fn,'w');
fprintf(fid,'SDDS1\n');
fprintf(fid,'&column name=t, units=s, type=double,  &end\n');
fprintf(fid,'&column name=Wz, units=V/C, type=double,  &end\n');	
fprintf(fid,'&data mode=ascii, &end\n');
fprintf(fid,'! page number 1\n');
fprintf(fid,'               %5.0f\n',length(pdata(:,1)));
for j = 1:length(pdata(:,1))
  fprintf(fid,'%12.8e  %12.8e\n',t(j),Wz(j));
end
fclose(fid);
