clc;clear
close all

fnum = 3;
for ss = 1:length(fnum)
    scanfold = 'BC1COL/col_';
    foldname = ['elegant/scan/',scanfold,num2str(fnum(ss))];
    ssnum = 51;
    rep_figure = 1;
    Q = 189.4225e-12;
    % plotname = {'BC2END','L3END','DL2END','UNDBEG'};
%     plotname = {'DL2END','UNDBEG','XLBEG'};
    plotname = {'XLBEG'};
    
    figurefold = ['elegant/scan/figures/',scanfold,num2str(fnum(ss))];
    if ~exist(figurefold,'dir')
        mkdir(figurefold);
    end
    
    for ii = 1:length(ssnum)
        subdatafold = [foldname,'/s',num2str(ssnum(ii))];
        subfigurefold = [figurefold,'/s',num2str(ssnum(ii))];
        if ~exist(subfigurefold,'dir')
            mkdir(subfigurefold);
        end
        for jj = 1:length(plotname)
            figurename = [subfigurefold,'/',plotname{jj},'.png'];
            if rep_figure || ~exist(figurename,'file')
                [X,Y,Z,I] = sddsplottp(subdatafold,plotname{jj},Q);
                figure(1)
                set(gcf,'position',[100,100,800,800])
                subplot(2,1,1)
                imagesc(X*1e15,Y/1000,Z)
                axis xy
                colormap(jetvar)
                xlabel('time (fs)')
                ylabel('E (GeV)')
                title(plotname{jj})
                set(gca,'xlim',[X(1)*1e15,X(end)*1e15]);
                enhance_plot('times',16,2,8)
                legend off
                yyaxis right
                plot(X*1e15,I,'r')
                ylabel('Current (kA)')
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
%                 print('-dpng','-r0',figurename)
            end
        end
    end
end