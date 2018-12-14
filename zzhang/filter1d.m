function y2 = filter1d(y,t,resl_t)
fft_y = fftshift(fft(y));
dt = t(2) - t(1);
DT = max(t) - min(t);
df = 1./DT;
fmax = 1./dt;
f = linspace(-fmax/2,fmax/2,length(y));

resl_f = 1./resl_t;

fft_y2 = fft_y.*exp(-f.^4/2/resl_f.^4);

y2 = ifft(ifftshift(fft_y2));

% figure;
% semilogy(f,abs(fft_y))
% hold on
% semilogy(f,abs(fft_y2))