function [fwhm,tc,tp,peak_power] = single_spike_analysis(t,power)
peak_power = max(power);
tp = mean(t(power == peak_power));
t_interp = linspace(min(t),max(t),length(t)*10);
power_interp = interp1(t,power,t_interp);
t_cut = t_interp(power_interp>0.5*peak_power);
fwhm = abs(t_cut(end) - t_cut(1));
t_range = fwhm/2.35*3;
power_cut2 = power_interp((t_interp>(tp-t_range))&(t_interp<(tp+t_range)));
t_cut2 = t_interp((t_interp>(tp-t_range))&(t_interp<(tp+t_range)));
tc = sum(t_cut2.*power_cut2)./sum(power_cut2);