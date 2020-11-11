close all
clc
clear all


axisticksize = 40;
origin = [200,60];
markersize = 20;
LineWidth = 2;
dbdown=40
map='default';
cut_pec=-50000%-0.1011%-0.42;
eks_pec=2;
below_pec=0;
noEn=1;
noport=1;
noEobs=1;
threshold=1000;




% noEobs=0;
% noport=0;
% noEn=0;
% nev=1; 
% name='pillbox_4000';
% freqstr='_freq_1145000000.0000000';

% 
% noEobs=0;
% noport=0;
% noEn=0;
% nev=1; 
% name='pillbox_50K';
% freqstr='_freq_1145000000.0000000';



% noport=0;
% name='cavity_rec_5K_feko';
% freqstr='_freq_1502400000.0000000';



noEobs=0;
noEn=0;
noport=0;
nev=20;
name='cavity_rec_17K_feko';
% freqstr='_freq_1510000000.0000000';
% freqstr='_freq_2327400000.0000000';
% freqstr='_freq_2316700000.0000000';
freqstr='_freq_2321000000.0000000';
% freqstr='_freq_7495000000.0000000';
% freqstr='_freq_2141300000.0000000';



% 
% noport=1;
% nev=100;   
% name='cavity_5cell_30K_feko';
% freqstr='_freq_937730000.00000000';
% freqstr='_freq_649960000.00000000';
% freqstr='_freq_648160000.00000000';
% freqstr='_freq_643620000.00000000';
% freqstr='_freq_638132000.00000000';
% freqstr='_freq_633770000.00000000';
% 



% name='pillbox_4000';

% 
% noport=0;
% name='cavity_wakefield_4K_feko';
% freqstr='_freq_1136100000.0000000.out';


fid2=fopen([name,'_node.inp'],'r');
fid3=fopen([name,'_elem.inp'],'r');


nnod2=fscanf(fid2,'%d',1);
npat=fscanf(fid3,'%d',1);

t1=cputime;
xyz2 = zeros(nnod2,3);
for i=1:nnod2
    tmp=fscanf(fid2,'%d',1);
    temp(1:3)=(fscanf(fid2,'%f',3))';
    xyz2(i,1)=temp(1);
    xyz2(i,2)=temp(2);
    xyz2(i,3)=temp(3);
end

t2=cputime;
display('surface node read time:')
t2-t1
t1=cputime;
ii=1;
eks=eks_pec;
cut=cut_pec;
below=below_pec;
ct2=zeros(npat,3);
if (below==1)
    for i=1:npat
        tmp=(fscanf(fid3,'%d',1))';
        ct2(ii,1:3)=(fscanf(fid3,'%d',3))';
%         if (xyz2(ct2(ii,1),eks)<=cut & xyz2(ct2(ii,2),eks)<=cut & xyz2(ct2(ii,3),eks)<=cut)
            ii=ii+1;
            %if any one of the nodes are outside the cut dont do anything
%         end
    end
else
    for i=1:npat
        tmp=(fscanf(fid3,'%d',1))';
        ct2(ii,1:3)=(fscanf(fid3,'%d',3))';
        %if (xyz2(ct2(ii,1),eks)>=cut & xyz2(ct2(ii,2),eks)>=cut & xyz2(ct2(ii,3),eks)>=cut ...
        %        & xyz2(ct2(ii,1),eks)<=5.999 & xyz2(ct2(ii,2),eks)<=5.999 & xyz2(ct2(ii,3),eks)<=5.999)
%         if (xyz2(ct2(ii,1),eks)>=cut & xyz2(ct2(ii,2),eks)>=cut & xyz2(ct2(ii,3),eks)>=cut)

            ii=ii+1;
            %if any one of the nodes are outside the cut dont do anything
%         end
    end
end
count_pec=ii-1
t2=cputime;
display('surface element read and processing time:')
t2-t1



ith=0;
for i=1:nev
if(noport==1)       
    fid1 = fopen (['code_outputs/EigVec_',num2str(i),freqstr,'.out'],'r');
else
    fid1 = fopen (['code_outputs/EigVecJ_',num2str(i),freqstr,'.out'],'r');
end
CV=zeros(nnod2,1);
% color for each node
CV(1:nnod2)=fscanf(fid1,'%f',nnod2)';

maxabs = max((abs(CV)));
CV = CV./maxabs;


if(norm(CV,1)>threshold)
    disp(['ith: ',num2str(i)])
    disp(['2-norm: ',num2str(norm(CV,2))])
    disp(['1-norm: ',num2str(norm(CV,1))])
    disp(['inf-norm: ',num2str(norm(CV,inf))])
    
    CV =20*log10(abs(CV));   
   
    X=xyz2(:,1);
    Y=xyz2(:,2);
    Z=xyz2(:,3);
    
    ith=ith+1;
    figure(ith);
    colormap(map);


    clf
    a = trisurf(ct2(1:count_pec,1:3),X,Y,Z,CV); %pause;

    c = colorbar('location','EastOutside');
    caxis([-dbdown 0]);
    axis equal;
    shading interp;
    axis off;
    grid off;
    view(-120,30);
    title(['$\mathbf{J}$ ',num2str(i),'th'],'interpreter','Latex')
    fname=sprintf('%s',name,'_currJ_',num2str(i),'th',freqstr,'.jpg');
    set(gcf,'Position',[origin,1000,700]);
    gca = get(gcf,'CurrentAxes');
    set( gca, 'FontName','Times New Roman','fontsize',axisticksize);    
    print ('-djpeg','-r300',fname)
    
    if(noport==0) 
   
        fid1 = fopen (['code_outputs/EigVecM_',num2str(i),freqstr,'.out'],'r');
        X=xyz2(:,1);
        Y=xyz2(:,2);
        Z=xyz2(:,3);

        CV=zeros(nnod2,1);
        % color for each node
        CV(1:nnod2)=fscanf(fid1,'%f',nnod2)';

        ith=ith+1;
        figure(ith);
        colormap(map);
        maxabs1 = max((abs(CV)));
        CV = CV./maxabs1;
        CV =20*log10(abs(CV));

        clf
        a = trisurf(ct2(1:count_pec,1:3),X,Y,Z,CV); %pause;

        c = colorbar('location','EastOutside');
        caxis([-dbdown 0]);
        axis equal;
        shading interp;
        axis off;
        grid off;
        view(-120,30);
        title(['$\mathbf{M}/\eta_0$ ',num2str(i),'th'],'interpreter','Latex')
        set(gcf,'Position',[origin,1000,700]);
        gca = get(gcf,'CurrentAxes');
        set( gca, 'FontName','Times New Roman','fontsize',axisticksize);        
        fname=sprintf('%s',name,'_currM_',num2str(i),'th',freqstr,'.jpg');
        print ('-djpeg','-r300',fname)


        fclose('all');

        axis equal
        shading interp
    end
    
    
     if(noEn==0) 
   
        fid1 = fopen (['code_outputs/EigEn_',num2str(i),freqstr,'.out'],'r');
        X=xyz2(:,1);
        Y=xyz2(:,2);
        Z=xyz2(:,3);

        CV=zeros(nnod2,1);
        % color for each node
        CV(1:nnod2)=fscanf(fid1,'%f',nnod2)';

        ith=ith+1;
        figure(ith);
        colormap(map);
        maxabs2 = max((abs(CV)));
        CV = CV./maxabs2;
        CV =20*log10(abs(CV));

        clf
        a = trisurf(ct2(1:count_pec,1:3),X,Y,Z,CV); %pause;

        c = colorbar('location','EastOutside');
        caxis([-dbdown 0]);
        axis equal;
        shading interp;
        axis off;
        grid off;
        view(-120,30);    
        title(['$|\mathbf{E}_n|$ ',num2str(i),'th'],'interpreter','Latex')
        set(gcf,'Position',[origin,1000,700]);
        gca = get(gcf,'CurrentAxes');
        set( gca, 'FontName','Times New Roman','fontsize',axisticksize);        

        fname=sprintf('%s',name,'_currM_',num2str(i),'th',freqstr,'.jpg');
        print ('-djpeg','-r300',fname)


        fclose('all');

        axis equal
        shading interp
     end   
    
    
     if(noEobs==0) 
   
        E_GF = load (['code_outputs/Eobs_',num2str(i),freqstr,'.out'],'r');
       
        ith=ith+1;
        figure(ith);
        plot(E_GF(:,3),E_GF(:,6),'r','LineWidth',LineWidth)
        gca = get(gcf,'CurrentAxes');
        set( gca, 'FontName','Times New Roman','fontsize',axisticksize);
        set(gcf,'Position',[origin,1000,700]);
        title(['$|\mathbf{E}|$ ',num2str(i),'th'],'interpreter','Latex')

     end 
     
     
     
end

end

