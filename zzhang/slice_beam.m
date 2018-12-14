function slicepara = slice_beam(x,xp,y,yp,t,p,Nbin,Q)
tmin = min(t);
tmax = max(t);
dt = (tmax - tmin)/Nbin;
tbin = (tmin-dt/2):dt:(tmax+dt/2);
q = Q/length(t);
Nbin = length(tbin);
I = zeros(1,Nbin-1);
avg_gamma = zeros(1,Nbin-1);
std_gamma = zeros(1,Nbin-1);
emitx = zeros(1,Nbin-1);
emity = zeros(1,Nbin-1);
betax = zeros(1,Nbin-1);
betay = zeros(1,Nbin-1);
alphax = zeros(1,Nbin-1);
alphay = zeros(1,Nbin-1);
sigx = zeros(1,Nbin-1);
sigy = zeros(1,Nbin-1);
xc = zeros(1,Nbin-1);
xpc = zeros(1,Nbin-1);
yc = zeros(1,Nbin-1);
ypc = zeros(1,Nbin-1);

for ii = 1:(Nbin-1)
    x_bin = x((t>tbin(ii))&(t<=tbin(ii+1)));
    xp_bin = xp((t>tbin(ii))&(t<=tbin(ii+1)));
    y_bin = y((t>tbin(ii))&(t<=tbin(ii+1)));
    yp_bin = yp((t>tbin(ii))&(t<=tbin(ii+1)));
    t_bin = t((t>tbin(ii))&(t<=tbin(ii+1)));
    p_bin = p((t>tbin(ii))&(t<=tbin(ii+1)));
    
    if ~isempty(t_bin)
        I(ii) = q*length(t_bin)/dt;
        avg_gamma(ii) = mean(p_bin);
        std_gamma(ii) = std(p_bin);
        xc(ii) = mean(x_bin);
        xpc(ii) = mean(xp_bin);
        yc(ii) = mean(y_bin);
        ypc(ii) = mean(yp_bin);
        x_bin = x_bin - mean(x_bin);
        xp_bin = xp_bin - mean(xp_bin);
        y_bin = y_bin - mean(y_bin);
        yp_bin = yp_bin - mean(yp_bin);
        emitx(ii) = sqrt(mean(x_bin.^2).*mean(xp_bin.^2)-mean(x_bin.*xp_bin).^2).*avg_gamma(ii);
        emity(ii) = sqrt(mean(y_bin.^2).*mean(yp_bin.^2)-mean(y_bin.*yp_bin).^2).*avg_gamma(ii);
        betax(ii) = mean(x_bin.*x_bin).*avg_gamma(ii)./emitx(ii);
        betay(ii) = mean(y_bin.*y_bin).*avg_gamma(ii)./emity(ii);
        alphax(ii) = -mean(x_bin.*xp_bin).*avg_gamma(ii)./emitx(ii);
        alphay(ii) = -mean(y_bin.*yp_bin).*avg_gamma(ii)./emity(ii);
        sigx(ii) = std(x_bin);
        sigy(ii) = std(y_bin);
    else
        I(ii) = 0;
        avg_gamma(ii) = 0;
        std_gamma(ii) = 0;
        emitx(ii) = 0;
        emity(ii) = 0;
        betax(ii) = 0;
        betay(ii) = 0;
        alphax(ii) = 0;
        alphay(ii) = 0;
        sigx(ii) = 0;
        sigy(ii) = 0;
    end
end

tbin2 = (tbin(1:(Nbin-1))+tbin(2:Nbin))/2;
slicepara.tbin = tbin2;
slicepara.current = I;
slicepara.avg_gamma = avg_gamma;
slicepara.std_gamma = std_gamma;
slicepara.emitx = emitx;
slicepara.emity = emity;
slicepara.betax = betax;
slicepara.betay = betay;
slicepara.alphax = alphax;
slicepara.alphay = alphay;
slicepara.sigx = sigx;
slicepara.sigy = sigy;
slicepara.xc = xc;
slicepara.xpc = xpc;
slicepara.yc = yc;
slicepara.ypc = ypc;
