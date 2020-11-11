clear all
clc
close all
E_GF = load('code_outputs/Eobs_1_freq_1145000000.0000000.out');
axisticksize = 40;
origin = [200,60];
markersize = 20;
LineWidth = 2;


figure(1)
plot(E_GF(:,3),E_GF(:,6),'r','LineWidth',LineWidth)
gca = get(gcf,'CurrentAxes');
set( gca, 'FontName','Times New Roman','fontsize',axisticksize);

set(gcf,'Position',[origin,1000,700]);