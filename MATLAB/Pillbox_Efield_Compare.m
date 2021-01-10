clear all
clc
close all
E_GF = load('code_outputs/fort.667');
E_LOC = load('code_outputs/fort.666');
E_REF = besselj(0,E_GF(:,1)*2.4048/0.1);
axisticksize = 40;
origin = [200,60];
markersize = 20;
LineWidth = 2;


figure(1)
plot(E_GF(:,1),E_GF(:,4)/max(E_GF(:,4)),'r','LineWidth',LineWidth)
hold on
plot(E_LOC(:,1),E_LOC(:,4)/max(E_GF(:,4)),'b','LineWidth',LineWidth)
hold on
plot(E_GF(:,1),E_REF/max(E_REF),'k','LineWidth',LineWidth)
hold on
legend('GF','Local','Reference')
gca = get(gcf,'CurrentAxes');
set( gca, 'FontName','Times New Roman','fontsize',axisticksize);

set(gcf,'Position',[origin,1000,700]);