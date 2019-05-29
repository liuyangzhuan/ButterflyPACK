close all;
clc;
clear all;
fontsize = 24;
axisticksize = 24;
origin = [200,60];
markersize = 20;
LineWidth = 3;


% rcs_new = load('bistaticH.out');
ith=5;
filename=['VV_bistatic_EigVec_',num2str(ith),'_freq_10600000000.000000.txt'];
rcs_new = load(filename);
plot(rcs_new(:,1),rcs_new(:,2),'k','LineWidth',LineWidth)
hold on


xlim([0,180])
set(gca,'xtick',[0,30,60,90,120,150,180]);
set( gca, 'FontName','Times New Roman','fontsize',axisticksize);
ylabel('Far-field Pattern [dB]')
xlabel('\it{\theta} \rm{[deg]}')
origin = [200,60];

grid on
gca = get(gcf,'CurrentAxes');
set(gca,'Xcolor',[0.8 0.8 0.8]); 
set(gca,'Ycolor',[0.8 0.8 0.8]); 
Caxes = copyobj(gca,gcf);
set(Caxes, 'color', 'none', 'xcolor', 'k', 'xgrid', 'off', 'ycolor','k', 'ygrid','off'); 
set(gcf,'Position',[origin,1000,700]);


fig = gcf;
style = hgexport('factorystyle');
style.Bounds = 'tight';
figname = ['far pattern_mode_',num2str(ith)];   
saveas(fig,figname,'meta')

