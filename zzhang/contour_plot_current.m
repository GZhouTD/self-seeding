function [X,Y,Z,dx,dy,I] = contour_plot_current(x,y,Nx,Ny,fig,Q,trange,prange)

%	function [X,Y,Z,dx,dy] = contour_plot(x,y,Nx,Ny,fig);
%
%	Function to generate 2D color image from 2D particle ditribution.
%
%	INPUTS:		x:		One of the two coordinates in the array of particles
%	     		y:		The other of the two coordinates in the array of particles
%				Nx:		The number of hor. bins
%				Ny:		The number of ver. bins
%				fig:	[Optional,DEF=1] if fig==1, get plot, else no plot,
%						just output arguments
%	OUTPUTS:	X:		The array of hor. bin-centers
%				Y:		The array of ver. bin-centers
%				Z:		The height of each point (pixel strength or color)
%				dx:		The calibration of x-input-units per bin (e.g., mm/bin)
%				dy:		The calibration of y-input-units per bin (e.g., mm/bin)

%==============================================================================

if ~exist('fig') %#ok<EXIST>
  fig = 1;
end

if ~exist('trange') %#ok<EXIST>
    xmin = min(x);
    xmax = max(x);
else
    xmin = trange(1);
    xmax = trange(2);
end

if ~exist('prange') %#ok<EXIST>
    ymin = min(y);
    ymax = max(y);
else
    ymin = prange(1);
    ymax = prange(2);
end

dx   = (xmax - xmin)/Nx;
dy   = (ymax - ymin)/Ny;
X    = ((xmin+dx/2):dx:(xmax-dx/2));
Y    = ((ymin+dy/2):dy:(ymax-dy/2));
Nx=length(X);
Ny=length(Y);

Z = zeros(Ny,Nx);
for j = 1:Ny
%     size(Y)
  i  = find(abs(y-Y(j))<dy/2);
  xi = x(i);
  yi = y(i);
  for k = 1:Nx
    ii = find(abs(xi-X(k))<dx/2);
    Z(j,k) = length(ii);
  end 
end
Np = length(x);
I = sum(Z,1)/dx*Q/Np/1000;

if fig>0
    if fig == 1
        figure
    else
        figure(fig)
    end
    subplot(2,1,1)
    imagesc(X*1e15,Y,Z)
    axis xy
    colormap(jetvar)
    xlabel('time (fs)')
    ylabel('\gamma')
    set(gca,'xlim',[X(1)*1e15,X(end)*1e15]);
    enhance_plot('times',16,2,8)
    legend off
    subplot(2,1,2)
    plot(X*1e15,I)
    xlabel('time (fs)')
    ylabel('Current (kA)')
    set(gca,'xlim',[X(1)*1e15,X(end)*1e15]);
    set(gca,'ylim',[0, max(I)*1.2]);
    enhance_plot('times',16,2,8)
    legend off
end