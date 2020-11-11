clear all;
close all;
clc;
% rcs_new = load('HH_bistatic.txt');
rcs_new = load('RCS_monostatic_ABCD_fix.txt');
plot(rcs_new(:,2),rcs_new(:,3),'k');
hold on

% rcs_new = load('RCS_monostatic_new_tfqmr.txt');
% % rcs_new = load('RCS_monostatic_new_tfqmr.txt');
% plot(rcs_new(:,2),rcs_new(:,3),'.r');
% hold on

rcs_new = load('RCS_monostatic_new_cori_intel17.txt');
plot(rcs_new(:,2),rcs_new(:,3),'r');
hold on

% rcs_new = load('RCS_monostatic_svd35.txt');
% plot(rcs_new(:,2),rcs_new(:,3),'b');
% hold on 


rcs_new = load('RCS_monostatic.txt'); 
plot(rcs_new(:,2),rcs_new(:,3),'g');
hold on


% 
% % % % rcs_new1 = load('VV_bistatic1.txt');
% plot(rcs_new(:,1),rcs_new(:,2));
hold on
% plot(rcs_new1(:,1),rcs_new1(:,2),'g');
% hold on

% rcs0 = load('RCS_monostatic_ref.txt');
% % rcs1 = load('../MLMDA_DIRECT_SOLVER_2D_EFIE_accuracy_reference/RCS_monostatic.txt');
% % rcs2 = load('../MLMDA_DIRECT_SOLVER_2D_EFIE_accuracy_reference/VV_bistatic2.txt');
% plot(rcs0(:,2),rcs0(:,3),'r');
% hold on
% plot(rcs1(:,2),rcs1(:,3),'k');
hold on
% plot(rcs2(:,1),rcs2(:,2),'k');
legend('ref nobplus','ref twolayer','new')
