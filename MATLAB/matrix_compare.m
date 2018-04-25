clear all;
close all;
clc;
axisticksize = 20;
origin = [200,60];
markersize = 8;
LineWidth = 2;



mm=1250;
nn=1250;
idx_s = nn*1146+1;
idx_e = nn*1150;
aa = load('fort.111');
bb = load('fort.112');
plot(aa(idx_s:idx_e))
hold on
plot(bb(idx_s:idx_e),'r')

mat_aa = zeros(mm,nn);
mat_bb = zeros(mm,nn);
for m=1:mm
    for n =1:nn
        mat_aa(m,n) = aa((m-1)*nn+n);
        mat_bb(m,n) = bb((m-1)*nn+n);        
    end
end

cmin = min(min(log10(mat_bb)));
cmax = max(max(log10(mat_bb)));

figure(3);
imagesc(log10(mat_aa));
axis equal
caxis([cmin,cmax]);
colorbar;
title('log(Z)')
gca = get(gcf,'CurrentAxes'); 
set( gca, 'FontName','Times New Roman','fontsize',axisticksize);
set(gcf,'Position',[origin,1000,700]);

figure(4);
imagesc(log10(mat_bb));
axis equal
caxis([cmin,cmax]);
colorbar;
title('log(Z)')
gca = get(gcf,'CurrentAxes'); 
set( gca, 'FontName','Times New Roman','fontsize',axisticksize);
set(gcf,'Position',[origin,1000,700]);


figure(5);
imagesc(log10(abs(mat_bb-mat_aa)./abs(mat_bb)));
axis equal
colorbar;
caxis([-10,0]);
gca = get(gcf,'CurrentAxes'); 
set( gca, 'FontName','Times New Roman','fontsize',axisticksize);
set(gcf,'Position',[origin,1000,700]);

norm(mat_bb-mat_aa)

