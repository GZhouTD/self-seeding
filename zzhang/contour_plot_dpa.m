function [X,Y,Z] = contour_plot_dpa(dpa,data,info,Nx,Ny)
cur = data.cur;
np = info.np;
dt = info.dt/3e8;
q = cur*dt/np;


x = dpa(:,5);
y = dpa(:,6);

xmin = min(x);
xmax = max(x);
ymin = min(y);
ymax = max(y);

dx   = (xmax - xmin)/Nx;
dy   = (ymax - ymin)/Ny;
X    = ((xmin+dx/2):dx:(xmax-dx/2));
Y    = ((ymin+dy/2):dy:(ymax-dy/2));
Nx=length(X);
Ny=length(Y);

Z = zeros(Ny,Nx);
queue = 1:length(x);
qq = q(ceil(queue/np));
for jj = 1:Ny
  xi  = x(abs(y-Y(jj))<dy/2);
  qi = qq(abs(y-Y(jj))<dy/2);
  for kk = 1:Nx
    qm = qi(abs(xi-X(kk))<dx/2);
    Z(jj,kk) = sum(qm);
  end 
end
figure
imagesc(X,Y,Z)
axis xy
colormap(jetvar)
xlabel('z (m)')
ylabel('\gamma')
enhance_plot('times',16,2,8)
