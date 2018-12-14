clc;clear
close all

load('scan_R56_65_high_peak_current.mat');

%%
[s1,s2] = size(output);
Imax = zeros(s1,s2);

R56_xleap1 = 100:10:400;
R56_xleap2 = 0:10:300;

for ii = 1:s1
    for jj = 1:s2
        [ii,jj]
%         figure(100)
%         set(gcf,'position',[100,100,800,800])
%         subplot(2,1,1)
%         imagesc(output{ii,jj}.X3*1e15,output{ii,jj}.Y3/1000,output{ii,jj}.Z3)
%         axis xy
%         colormap(jetvar)
%         xlabel('time (fs)')
%         ylabel('E (GeV)')
%         title(['R1 = ',num2str(R56_xleap1(ii)),'  R2 = ',num2str(R56_xleap2(jj))]);
%         enhance_plot('times',16,2,8)
%         legend off
%         subplot(2,1,2)
%         plot(output{ii,jj}.X3*1e15,output{ii,jj}.I3)
%         xlabel('time (fs)')
%         ylabel('Current (kA)')
%         set(gca,'xlim',[output{ii,jj}.X3(1)*1e15,output{ii,jj}.X3(end)*1e15]);
%         set(gca,'ylim',[0, max(output{ii,jj}.I3)*1.2]);
%         enhance_plot('times',16,2,8)
%         legend off
        
        I = output{ii,jj}.I3;
        t = output{ii,jj}.X3*1e15;
        
        % 65
        I((t<7)|(t>16)) = [];
        t((t<7)|(t>16)) = [];
        
        % 63
%         I((t<-3)|(t>7)) = [];
%         t((t<-3)|(t>7)) = [];
        % 62s
%         I((t<-8)|(t>3)) = [];
%         t((t<-8)|(t>3)) = [];
       % 51
%         I((t<-20.8)|(t>-10)) = [];
%         t((t<-20.8)|(t>-10)) = [];
        Imax(ii,jj) = max(I);
%         pause
    end
end

figure
imagesc(R56_xleap1,R56_xleap2,Imax')
axis xy
xlabel('R56 (\mum)')
ylabel('R56 (\mum)')
enhance_plot('times',16,2,8)
legend off
