clear all;
close all;
clc;
% rcs_new = load('bistaticH.out');
rcs_new = load('HH_bistatic.txt');
plot(rcs_new(:,1),rcs_new(:,2));
hold on
% plot(rcs_new1(:,1),rcs_new1(:,2),'g');
% hold on




% % rcs0 = load('../MLMDA_DIRECT_SOLVER_2D_EFIE_accuracy_reference/VV_bistatic.txt');
% % % rcs1 = load('../MLMDA_DIRECT_SOLVER_2D_EFIE_accuracy_reference/VV_bistatic1.txt');
% % % rcs2 = load('../MLMDA_DIRECT_SOLVER_2D_EFIE_accuracy_reference/VV_bistatic2.txt');
% % plot(rcs0(:,1),rcs0(:,2),'r');



% hold on
% plot(rcs1(:,1),rcs1(:,2),'k');
% hold on
% plot(rcs2(:,1),rcs2(:,2),'k');
legend('new','reference')
