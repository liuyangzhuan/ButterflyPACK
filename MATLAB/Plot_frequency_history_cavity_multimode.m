close all
clc
clear all


axisticksize = 28;
origin = [200,60];
markersize = 10;
LineWidth = 3;


name='cavity_rec_17K_feko';


fid1=fopen(['code_outputs/',name,'_Nmodes.txt'],'r');
nmode=fscanf(fid1,'%d',1);
fclose(fid1)

lineColors = line_colors(nmode);



fid1=fopen(['code_outputs/',name,'_freq_history.txt'],'r');
nfreq=fscanf(fid1,'%d',1);
freqs=zeros(nfreq,1);
for i=1:nfreq
    tmp=fscanf(fid1,'%f',1);
    freqs(i)=round(tmp/1e5);
end    
fclose(fid1)
freqs=sort(freqs);

eigval0 = ones(nfreq,1);

nthmode = 1:nmode;
for nn=1:nmode
    hist = load(['code_outputs/',name,'_EigVals_',num2str(nn),'.out'])
    [~,idx] = sort(hist(:,1));
    hist_sort=hist(idx,:)
    hist_sort(:,1)=round(hist_sort(:,1)/1e5);
    nthmode(nn)= hist_sort(find(hist_sort(:,2)==min(hist_sort(:,2))),1);
end
[~,idxmode] = sort(nthmode);

leg={};
for tt=1:nmode
    leg{tt}=['mode ',num2str(tt)];
    nn=idxmode(tt)
    hist = load(['code_outputs/',name,'_EigVals_',num2str(nn),'.out'])
    [~,idx] = sort(hist(:,1));
    hist_sort=hist(idx,:)
    hist_sort(:,1)=round(hist_sort(:,1)/1e5);
    eigval = eigval0;
    
    for ii=1:length(hist_sort)
       eigval(find(freqs(:,1)==hist_sort(ii,1)))= hist_sort(ii,2);
    end
    

%     figure(nn)
    hold on
%     plot(freqs/1e4,eigval,'o-','MarkerSize',markersize,'MarkerFaceColor',[lineColors(nn,:)],'LineWidth',LineWidth, 'Color',[lineColors(nn,:)])
    plot(hist_sort(:,1)/1e4,hist_sort(:,2),'o-','MarkerSize',markersize,'MarkerFaceColor',[lineColors(nn,:)],'LineWidth',LineWidth, 'Color',[lineColors(nn,:)])
    set(gca,'yscale','log') 
end


format short
disp('The resonant frequencies (GHz):')
nthmode(idxmode)/1e4


gca = get(gcf,'CurrentAxes');
set( gca, 'FontName','Times New Roman','fontsize',axisticksize);
set(gcf,'Position',[origin,1000,700]);
ylabel('Min Eigenvalue','interpreter','Latex')    
xlabel('Frequency (GHz)','interpreter','Latex')  
legend(leg,'Location','northeast','color','none')
fig = gcf;

str = 'history_cavity.pdf';

saveas(fig,str,'epsc')
print(str,'-dpdf','-bestfit')






