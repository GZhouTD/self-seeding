clc;clear
close all

foldname = 'elegant/scan/BC1COL/col_3';
ssnum = 2;
charge = 189.4225e-12;

Q = charge;
plotname = 'UNDBEG';

for ii = 1:length(ssnum)
    sddsfile = [foldname,'/s',num2str(ssnum(ii)),'/',plotname,'.out'];
    plainfile = [foldname,'/s',num2str(ssnum(ii)),'/',plotname,'.pla'];
    covert_command = ['sdds2plaindata ',sddsfile,' ',plainfile,' -col=t -col=p -sep=''  '''];
    if ~exist(plainfile,'file')
        system(covert_command);
    end
    fid = fopen(plainfile);
    np = fscanf(fid,'%d',1);
    pdata = fscanf(fid,'%f',np*2);
    fclose(fid);
    t = pdata(1:2:np*2);
    p = pdata(2:2:np*2);
    t = t - mean(t);
    [X,Y,Z,~,~,I] = contour_plot_current(t,p*0.511,300,300,0,Q);
    figure(2)
    set(gcf,'position',[100,100,800,800])
    subplot(2,1,1)
    imagesc(X*1e15,Y/1000,Z)
    axis xy
    colormap(jetvar)
    xlabel('time (fs)')
    ylabel('E (GeV)')
    title([plotname])
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
    [tbin,Eavg,Erms] = sliceEnergy(t,p,300);
    figure
    plot(tbin/1e-15,Erms*0.511)
end