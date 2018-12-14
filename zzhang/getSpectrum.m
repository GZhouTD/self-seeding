function [spectrum]=getSpectrum(amplitude,phase)
% calculates the spectrum out of a given intensity and phase

nz=size(amplitude,1);
nt=size(amplitude,2);

spectrum=zeros(nz,nt);


for i=1:nz
  signal=sqrt(amplitude(i,:)).*complex(cos(phase(i,:)),sin(phase(i,:)));
  spectrum(i,:)=abs(fftshift(fft(signal)));
  spectrum(i,:)=spectrum(i,:).*spectrum(i,:);
end
