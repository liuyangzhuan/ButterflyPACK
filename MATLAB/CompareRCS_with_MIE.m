close all;
clc;
clear all;
fontsize = 24;
axisticksize = 24;
origin = [200,60];
markersize = 10;
LineWidth = 3;



rcs_hh = load('../build/HH_bistatic.txt');
rcs_vv = load('../build/VV_bistatic.txt');
% E2 = 10.^(rcs_vv(:,2)/10)+10.^(rcs_hh(:,2)/10);


%% Calculate and plot the analytical results 
cc = 299792458.0e0;
lamda = 0.14;
f = cc/lamda;
omega = 2*pi*f;
k = omega/cc;
a = 1;
r = 100*lamda;
ka =k*a;
kr =k*r;
Nmax = 100;
theta = linspace(0,pi,1000);
RCS_cal_anl_hh = zeros(2,length(theta));
RCS_cal_anl_vv = zeros(2,length(theta));

an = zeros(Nmax,1);
bn = zeros(Nmax,1);
for nn=1:Nmax
    an(nn) = sphj(nn,ka)/sphh(nn,ka);
    bn(nn) = (ka*sphj(nn-1,ka)-nn*sphj(nn,ka))/(ka*sphh(nn-1,ka) - nn*sphh(nn,ka));
end
nn = 1:Nmax;
for idxtheta = 1:length(theta)
    phi = 0;
    Pkdtmp = legv1(Nmax,cos(theta(idxtheta)));
    Pkd = Pkdtmp(1,2:end);
    Pkd_hat = Pkdtmp(2,2:end);
    Pkd = Pkd/sin(theta(idxtheta));
    Pkd_hat = -Pkd_hat*sin(theta(idxtheta));
    Etheta(idxtheta) = j*exp(-j*ka)*cos(phi)/kr   *   sum((-1).^(nn).*(2*nn+1)./(nn.*(nn+1)).*(bn'.*Pkd_hat - an'.*Pkd));
    Ephi(idxtheta) = -j*exp(-j*ka)*sin(phi)/kr    *   sum((-1).^(nn).*(2*nn+1)./(nn.*(nn+1)).*(bn'.*Pkd - an'.*Pkd_hat));
    RCS_cal_anl_hh(:,idxtheta) = [Etheta(idxtheta); Ephi(idxtheta)];

    Etheta(idxtheta) = j*exp(-j*ka)*cos(phi)/kr   *   sum((-1).^(nn).*(2*nn+1)./(nn.*(nn+1)).*(an'.*Pkd_hat - bn'.*Pkd));
    Ephi(idxtheta) = -j*exp(-j*ka)*sin(phi)/kr    *   sum((-1).^(nn).*(2*nn+1)./(nn.*(nn+1)).*(an'.*Pkd - bn'.*Pkd_hat));
    RCS_cal_anl_vv(:,idxtheta) = [Etheta(idxtheta); Ephi(idxtheta)];

end 


subplot(2, 2, 1);
plot(theta*180/pi,10*log10(sum(abs(RCS_cal_anl_hh).^2)*4*pi*r.^2),'r','LineWidth',LineWidth)
hold on
plot(rcs_hh(:,1),10*log10((rcs_hh(:,2).^2+rcs_hh(:,3).^2)/4/pi),'k.','MarkerSize',markersize,'MarkerFaceColor','k')

set( gca, 'FontName','Times New Roman','fontsize',axisticksize);
legend({'Mie','ButterflyPACK'},'Color', 'none')
ylabel('bistatic RCS-hh (dB)')
xlabel('\phi')
xlim([0,180])
origin = [200,60];

grid on
% gca = get(gcf,'CurrentAxes');
% set(gca,'Xcolor',[0.8 0.8 0.8]); 
% set(gca,'Ycolor',[0.8 0.8 0.8]); 
% Caxes = copyobj(gca,gcf);
% set(Caxes, 'color', 'none', 'xcolor', 'k', 'xgrid', 'off', 'ycolor','k', 'ygrid','off'); 
set(gcf,'Position',[origin,1000,700]);
% % pause(2);
% fig = gcf;
% figname='RCS_BPACK_vs_Mie_mag_hh';
% saveas(fig,figname,'png')




subplot(2, 2, 2);
plot(theta*180/pi,rad2deg(unwrap(angle(RCS_cal_anl_hh(1,:)+RCS_cal_anl_hh(2,:)))),'r','LineWidth',LineWidth)
hold on
plot(rcs_hh(:,1),rad2deg(unwrap(angle(1i*rcs_hh(:,2) + rcs_hh(:,3)))),'k.','MarkerSize',markersize,'MarkerFaceColor','k')
gca = get(gcf,'CurrentAxes');
set( gca, 'FontName','Times New Roman','fontsize',axisticksize);
legend({'Mie','ButterflyPACK'},'Color', 'none')
ylabel('bistatic RCS-hh (degree)')
xlabel('\phi')
xlim([0,180])
origin = [200,60];

grid on
% gca = get(gcf,'CurrentAxes');
% set(gca,'Xcolor',[0.8 0.8 0.8]); 
% set(gca,'Ycolor',[0.8 0.8 0.8]); 
% Caxes = copyobj(gca,gcf);
% set(Caxes, 'color', 'none', 'xcolor', 'k', 'xgrid', 'off', 'ycolor','k', 'ygrid','off'); 
set(gcf,'Position',[origin,1000,700]);
% % pause(2);
% fig = gcf;
% figname='RCS_BPACK_vs_Mie_phase_hh';
% saveas(fig,figname,'png')



subplot(2, 2, 3);
plot(theta*180/pi,10*log10(sum(abs(RCS_cal_anl_vv).^2)*4*pi*r.^2),'r','LineWidth',LineWidth)
hold on
plot(rcs_vv(:,1),10*log10((rcs_vv(:,2).^2+rcs_vv(:,3).^2)/4/pi),'k.','MarkerSize',markersize,'MarkerFaceColor','k')

gca = get(gcf,'CurrentAxes');
set( gca, 'FontName','Times New Roman','fontsize',axisticksize);
legend({'Mie','ButterflyPACK'},'Color', 'none')
ylabel('bistatic RCS-vv (dB)')
xlabel('\phi')
xlim([0,180])
origin = [200,60];

grid on
% gca = get(gcf,'CurrentAxes');
% set(gca,'Xcolor',[0.8 0.8 0.8]); 
% set(gca,'Ycolor',[0.8 0.8 0.8]); 
% Caxes = copyobj(gca,gcf);
% set(Caxes, 'color', 'none', 'xcolor', 'k', 'xgrid', 'off', 'ycolor','k', 'ygrid','off'); 
set(gcf,'Position',[origin,1000,700]);
% % pause(2);
% fig = gcf;
% figname='RCS_BPACK_vs_Mie_mag_vv';
% saveas(fig,figname,'png')



subplot(2, 2, 4);
plot(theta*180/pi,rad2deg(unwrap(angle(RCS_cal_anl_vv(1,:)+RCS_cal_anl_vv(2,:)))),'r','LineWidth',LineWidth)
hold on
plot(rcs_vv(:,1),rad2deg(unwrap(angle(1i*rcs_vv(:,2) + rcs_vv(:,3))))-180,'k.','MarkerSize',markersize,'MarkerFaceColor','k')
gca = get(gcf,'CurrentAxes');
set( gca, 'FontName','Times New Roman','fontsize',axisticksize);
legend({'Mie','ButterflyPACK'},'Color', 'none')
ylabel('bistatic RCS-vv (degree)')
xlabel('\phi')
xlim([0,180])
origin = [200,60];

grid on
% gca = get(gcf,'CurrentAxes');
% set(gca,'Xcolor',[0.8 0.8 0.8]); 
% set(gca,'Ycolor',[0.8 0.8 0.8]); 
% Caxes = copyobj(gca,gcf);
% set(Caxes, 'color', 'none', 'xcolor', 'k', 'xgrid', 'off', 'ycolor','k', 'ygrid','off'); 
set(gcf,'Position',[origin,1000,700]);
% pause(2);
fig = gcf;
figname='RCS_BPACK_vs_Mie';
saveas(fig,figname,'png')