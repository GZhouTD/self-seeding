function y = smooth1d(y,ftime)
n = length(y);
for ss = 1:ftime
    for ii = 3:(n-2)
        y(ii) = 0.1*y(ii-2)+0.2*y(ii-1)+0.4*y(ii)+0.2*y(ii+1)+0.1*y(ii+2);
    end
end