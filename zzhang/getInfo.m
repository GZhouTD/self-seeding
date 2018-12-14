function getInfo
% display basic information about output data in memory

global output info

if (info.isloaded==0)
  fprintf('No output data into memory\n');
  return
end

fprintf('File Name: %s\n',info.file);
if (info.nt == 1)
    fprintf('Run: Steady-State\n');
elseif (info.isscan)
    fprintf('Run: Scan\n');
else
    fprintf('Run: Time-dependent\n');
end    
fprintf('Steps in z: %d\n',info.nz);
fprintf('Undulator ends at z = %f\n',output.z(info.nz));
if (info.nt > 1)
  fprintf('Steps in s: %d\n',info.nt);
  if (info.isscan ==0)
     fprintf('Bunch ends at s = %g\n',output.t(info.nt));
     fprintf('Slice spacing: %g\n',info.dt); 
  end
end
fprintf('Reference wavelength: %g\n',info.lambda);

fprintf('Output available for: \n');

if (size(output.power,1) > 1) fprintf(' power\n'); end
if (size(output.increment,1) > 1) fprintf(' increment\n'); end
if (size(output.signal,1) > 1) fprintf(' signal\n'); end
if (size(output.signalphase,1) > 1) fprintf(' signalphase\n'); end
if (size(output.radsize,1) > 1) fprintf(' radsize\n'); end
if (size(output.divergence,1) > 1) fprintf(' divergence\n'); end
if (size(output.energy,1) > 1) fprintf(' energy\n'); end
if (size(output.bunching,1) > 1) fprintf(' bunching\n'); end
if (size(output.xsize,1) > 1) fprintf(' xsize\n'); end
if (size(output.ysize,1) > 1) fprintf(' ysize\n'); end
if (size(output.xpos,1) > 1) fprintf(' xpos\n'); end
if (size(output.ypos,1) > 1) fprintf(' ypos\n'); end
if (size(output.energyspread,1) > 1) fprintf(' energyspread\n'); end
if (size(output.farfield,1) > 1) fprintf(' farfield\n'); end
if (size(output.spectrum,1) > 1) fprintf(' spectrum\n'); end
if (size(output.bandwidth,2) > 1) fprintf(' bandwidth\n'); end
if (size(output.error,1) > 1) fprintf(' error\n'); end
if (size(output.cur,2) > 1) fprintf(' current\n'); end


