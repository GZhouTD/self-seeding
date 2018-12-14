function [ output,info ] = readOutput( filename )
% filename='rt_seg_adsame_00.out';
%reads genesis output file into memory
% function takes pathname and filename to open file and parse the input
% processed data are stored into memory

% global info output

info=struct('isloaded',0,'file','','nz',0,'nt',0,'lambda',0.0,'dt',0.0,...
    'ng',0,'dg',0.0,'np',0,'npz',0,'npt',0,'nfz',0,'nft',0,'isscan',0);

output=struct('z',[0.],'t',[0.],'freq',[0.],'aw',[0.],'qf',[0.],...
    'power',[0.],'increment',[0.],'signal',[0.],'signalphase',[0.],...
    'radsize',[0.],'divergence',[0.],'energy',[0.],'bunching',[0.],...
    'xsize',[0.],'ysize',[0.],'error',[0.],'xpos',[0.],'ypos',[0.],...
    'energyspread',[0.],'farfield',[0.],'cur',[0.],'spectrum',[0.],...
    'bandwidth',[0.],'bunspec',[0.],'bun2',[0.],'pow2',[0.],'bun3',...
    [0.],'pow3',[0.],'bun4',[0.],'pow4',[0.],'bun5',[0.],'pow5',...
    [0.]);

fid=fopen(filename,'r');
if (fid < 1)
    err=-1;
    return;
end

tline=fgetl(fid);
while (ischar(tline)) 
  if (strfind(tline,'$end')>0)
      break;
  end    
  tline=fgetl(fid);  
end

if (~ischar(tline))  % end of file without finding '$end'
  err = -2;
  return;
end  

info.file=filename;

tline=fgetl(fid); % skip one line

% get the lout array
[lout,count] = fscanf(fid,' %d',19);
tline=fgetl(fid); % skip rest of line  
lsize=size(lout);
flag_main=lout(1:15);
flag_harm=lout(16:lsize(1));
record_column=sum(flag_main)+sum(flag_harm); % sorting of iput fields is done later
lout_modified=[flag_main;ones(length(flag_harm)*4,1)];

%header information, stored in the info struct
info.nz=fscanf(fid, '%d',1);   % number of integration steps
tline=fgetl(fid);
info.nt=fscanf(fid, '%d',1);    % numbers of records
tline=fgetl(fid);
info.lambda=fscanf(fid, '%e',1); % reference wavelength
tline=fgetl(fid);
info.dt=fscanf(fid, '%e',1);   % spacing in time-dependent simulation
tline=fgetl(fid);
info.ng=fscanf(fid, '%d',1);    % numbers of grid points in x or y
tline=fgetl(fid);
info.dg=fscanf(fid, '%e',1);   % grid spacing in x or y
tline=fgetl(fid);
info.np=fscanf(fid, '%d',1);   % number of particles per record
tline=fgetl(fid);
info.npz=fscanf(fid, '%d',1);   % number of particle records in z
tline=fgetl(fid);
info.npt=fscanf(fid, '%d',1);    % numbers of particle records in t
tline=fgetl(fid);
info.nfz=fscanf(fid, '%d',1);   % number of field records in z
tline=fgetl(fid); 
info.nft=fscanf(fid, '%d',1);    % numbers of field records in t
tline=fgetl(fid);

tline=fgetl(fid);    % skip one line
%read magnetic lattice
zglobal=fscanf(fid, '%e',[3 info.nz]);  % read 3 columns and nz rows
tline=fgetl(fid) ; % skip last of the line

data=zeros(record_column,info.nz,info.nt);
output.cur=[1:info.nt]*1.;

info.isscan=0;

for i =1:info.nt
  % skip 3 lines
  tline=fgetl(fid);
  tline=fgetl(fid);
  tline=fgetl(fid);
  tlocal=fscanf(fid,'%e',1);   % get current or scan parameter
  output.cur(i:i)=tlocal; 
  if (strfind(tline,'scan')>0)
    info.isscan=1;
  end  
  tline=fgetl(fid); % skip rest of line
  % skip 3 lines
  tline=fgetl(fid);
  tline=fgetl(fid);
  tline=fgetl(fid);
  record=fscanf(fid, '%e',[record_column info.nz]);  % read columns and nz rows
  data(:,:,i)=record; 
  tline=fgetl(fid); %skip rest of line
end
fclose(fid);

output.t=[1:info.nt]*info.dt; 
output.z=zglobal(1,:);
output.aw=zglobal(2,:);
output.qf=zglobal(3,:);
freq0=299792458/info.lambda;
if(info.isscan > 0)
   output.t=output.cur;
else  
  if (info.nt > 1)
    dfreq=299792458/(info.nt-1)/info.dt;
    output.freq=[0:info.nt-1]*dfreq-0.5*(info.nt-1.)*dfreq+freq0; % frequency
    output.freq=output.freq*6.62606957e-34/1.60217657e-19;
  end    
end

icount=0;
for i=1:length(lout_modified)
  if (lout_modified(i) ~= 0)
     icount=icount+1; 
     switch i
         case 1
           output.power=reshape(data(icount,:,:),info.nz,info.nt);
         case 2  
           output.increment=reshape(data(icount,:,:),info.nz,info.nt);
         case 3
           output.signal=reshape(data(icount,:,:),info.nz,info.nt);
         case 4
           output.signalphase=reshape(data(icount,:,:),info.nz,info.nt);
         case 5 
            output.radsize=reshape(data(icount,:,:),info.nz,info.nt);
         case 6
             output.divergence=reshape(data(icount,:,:),info.nz,info.nt);
         case 7
             output.energy=reshape(data(icount,:,:),info.nz,info.nt);
         case 8
             output.bunching=reshape(data(icount,:,:),info.nz,info.nt);
         case 9
             output.xsize=reshape(data(icount,:,:),info.nz,info.nt);
         case 10
             output.ysize=reshape(data(icount,:,:),info.nz,info.nt);
         case 11
             output.error=reshape(data(icount,:,:),info.nz,info.nt);
         case 12
             output.xpos=reshape(data(icount,:,:),info.nz,info.nt);
         case 13
             output.ypos=reshape(data(icount,:,:),info.nz,info.nt);
         case 14
             output.energyspread=reshape(data(icount,:,:),info.nz,info.nt);
         case 15
             output.farfield=reshape(data(icount,:,:),info.nz,info.nt);
%          case 16
%              output.bun2=reshape(data(icount,:,:),info.nz,info.nt);
%          case 17
%              output.pow2=reshape(data(icount,:,:),info.nz,info.nt);
%          case 18
%              
%          case 19
%              
%          case 20
%              output.bun3=reshape(data(icount,:,:),info.nz,info.nt);
%          case 21
%              output.pow3=reshape(data(icount,:,:),info.nz,info.nt);
%          case 22
%              
%          case 23
%              
%          case 24
%              output.bun4=reshape(data(icount,:,:),info.nz,info.nt);
%          case 25
%              output.pow4=reshape(data(icount,:,:),info.nz,info.nt);
%          case 26
%              
%          case 27
%              
%          case 28
%              output.bun5=reshape(data(icount,:,:),info.nz,info.nt);
%          case 29
%              output.pow5=reshape(data(icount,:,:),info.nz,info.nt);
%          case 30
%              
%          case 31
%              
     end
  end    
end

% calculate spectrum
if (info.isscan==0) && (info.nt > 1) && (size(output.signal,2) > 1) && (size(output.signalphase,2) >1 )
   output.spectrum=getSpectrum(output.signal,output.signalphase);
   output.bunspec=getSpectrum(output.bunching,output.signalphase);
   output.bandwidth=[1:info.nz]*1.;
   for i=1:info.nz
      w1=sum(output.spectrum(i,:));
      w2=sum(output.freq.*output.spectrum(i,:));
      w3=sum(output.freq.*output.freq.*output.spectrum(i,:));
      output.bandwidth(i)=sqrt(w3/w1-w2*w2/w1/w1);
   end
   output.bandwidth=output.bandwidth;
end

info.isloaded=1;

clear data;




