clear all;
clc;
close all;

axisticksize = 24;
origin = [200,60];
markersize = 2;
LineWidth = 2;


aa = load('IterHistory.txt_r12_ButtLevel3_LS5');
plot(aa(:,1),log10(aa(:,2)),'b','LineWidth',LineWidth);
hold on
aa = load('IterHistory.txt_r15_ButtLevel3_LS5');
plot(aa(:,1),log10(aa(:,2)),'r','LineWidth',LineWidth);
hold on
aa = load('IterHistory.txt_r20_ButtLevel3_LS5');
plot(aa(:,1),log10(aa(:,2)),'k','LineWidth',LineWidth);
hold on


% aa = load('IterHistory.txt_r12_ButtLevel3_LS7');
% plot(aa(:,1),log10(aa(:,2)),'b','LineWidth',LineWidth);
% hold on
% aa = load('IterHistory.txt_r15_ButtLevel3_LS7');
% plot(aa(:,1),log10(aa(:,2)),'r','LineWidth',LineWidth);
% hold on
% aa = load('IterHistory.txt_r20_ButtLevel3_LS7');
% plot(aa(:,1),log10(aa(:,2)),'k','LineWidth',LineWidth);
% hold on


% aa = load('IterHistory.txt_r12_ButtLevel3_LS9');
% plot(aa(:,1),log10(aa(:,2)),'b','LineWidth',LineWidth);
% hold on
% aa = load('IterHistory.txt_r15_ButtLevel3_LS9');
% plot(aa(:,1),log10(aa(:,2)),'r','LineWidth',LineWidth);
% hold on
% aa = load('IterHistory.txt_r20_ButtLevel3_LS9');
% plot(aa(:,1),log10(aa(:,2)),'k','LineWidth',LineWidth);
% hold on

gca = get(gcf,'CurrentAxes'); 
set( gca, 'FontName','Times New Roman','fontsize',axisticksize);

legend('rank=12','rank=15','rank=20',2)
xlabel('Iteration Count')
ylabel('Error')
ylim([-7,0])

gca = get(gcf,'CurrentAxes');
set(gca,'Xcolor',[0.8 0.8 0.8]); 
set(gca,'Ycolor',[0.8 0.8 0.8]); 
Caxes = copyobj(gca,gcf);
set(Caxes, 'color', 'none', 'xcolor', 'k', 'xgrid', 'off', 'ycolor','k', 'ygrid','off'); 
set(gcf,'Position',[origin,1000,700]);