close all;
clear all;
clc;
format shortE
FSize=50;
Fweight='normal';
load './seismic.map'
twoconst=[''];


lineColors = line_colors(6);
% lineColors(1,:)=[7.2000e-01   8.0000e-01   4.0000e-01];

% h2s = [0.01 0.005];
% yobs = [0.25 0.5-max(h2s)*10];
% zobs = [0.25 0.5-max(h2s)*10];
% f=5;
% vs=1;
% ivelo=1;ex=1;shape=-1;
% smax=2;
% tol=1e-5;  % tolerence in weno_bf
% h1=h2s(1);  %0.015625; %0.015625; % grid size in weno_bf
% h3=h2s(1);  % grid size in 9-point fdfd




% h2s = [0.005 0.0025];
% yobs = [0.25 0.5-max(h2s)*10];
% zobs = [0.25 0.5-max(h2s)*10];
% f=20;
% vs=1;
% ivelo=1;ex=1;shape=-1;
% smax=2;
% tol=1e-6;  % tolerence in weno_bf
% h1=h2s(1);  %0.015625; %0.015625; % grid size in weno_bf
% h3=h2s(1);  % grid size in 9-point fdfd







% h2s = [0.01 0.005];
% yobs = [0.25 0.5-max(h2s)*10];
% zobs = [0.25 0.5-max(h2s)*10];
% f=4;
% vs=1;
% ivelo=9;ex=2;shape=-1;
% smax=3.3333;
% tol=1e-5;  % tolerence in weno_bf
% h1=h2s(1);  %0.015625; %0.015625; % grid size in weno_bf
% h3=h2s(1);  % grid size in 9-point fdfd
% 

% 
% h2s = [0.004 0.002];
% % yobs = [0.6 0.52-max(h2s)*10];
% % zobs = [0.6 0.52-max(h2s)*10];
% yobs = [0.6 max(h2s)*10];
% zobs = [0.6 max(h2s)*10];
% 
% f=16;
% vs=1;
% ivelo=9;ex=2;shape=-1;
% smax=3.42466;
% tol=1e-5;  % tolerence in weno_bf
% h1=h2s(1);  %0.015625; %0.015625; % grid size in weno_bf
% h3=h2s(1);  % grid size in 9-point fdfd
% 
% 

% 
% h2s = [0.02];
% % yobs = [0.6 0.52-max(h2s)*10];
% % zobs = [0.6 0.52-max(h2s)*10];
% yobs = [0.6 max(h2s)*10];
% zobs = [0.6 max(h2s)*10];
% 
% f=4;
% vs=1;
% ivelo=9;shape=4;
% smax=3.42466;
% tol=1e-5;  % tolerence in weno_bf
% h1=h2s(1);  %0.015625; %0.015625; % grid size in weno_bf
% % twoconst=['_ivelo9_const'];
% twoconst=[''];
% 




% 
h2s = [0.02];
% yobs = [0.6 0.52-max(h2s)*10];
% zobs = [0.6 0.52-max(h2s)*10];
yobs = [0.6 max(h2s)*10];
zobs = [0.6 max(h2s)*10];

f=4;
vs=1;
ivelo=9;shape=4;
smax=8.33333;
tol=1e-4;  % tolerence in weno_bf
h1=h2s(1);  %0.015625; %0.015625; % grid size in weno_bf
h3=h2s(1);  % grid size in 9-point fdfd
L=0.8;
center=[0.4,0.4,0.4];


% % % 
% 
% h2s = [0.01];
% % yobs = [0.6 0.52-max(h2s)*10];
% % zobs = [0.6 0.52-max(h2s)*10];
% yobs = [0.6 max(h2s)*10];
% zobs = [0.6 max(h2s)*10];
% 
% f=8;
% vs=1;
% ivelo=11;shape=4;
% smax=3;
% tol=1e-5;  % tolerence in weno_bf
% h1=h2s(1);  %0.015625; %0.015625; % grid size in weno_bf
% % twoconst=['_ivelo9_const'];
% twoconst=[''];
% L=0.8;
% center=[0.4,0.4,0.4];





%% Draw the bounding box for the inhomogenous domain
xMin = center(1)-L/2;  xMax = center(1)+L/2;
yMin = center(2)-L/2;  yMax = center(2)+L/2;
zMin = center(3)-L/2;  zMax = center(3)+L/2;
% Corner (vertex) coordinates:
vertices = [
xMin  yMin  zMin
xMax  yMin  zMin
xMax  yMax  zMin
xMin  yMax  zMin
xMin  yMin  zMax
xMax  yMin  zMax
xMax  yMax  zMax
xMin  yMax  zMax
];
% Faces (each row = one quadrilateral given by 4 vertices):
faces = [
1 2 3 4   % bottom
5 6 7 8   % top
1 2 6 5   % side
2 3 7 6
3 4 8 7
4 1 5 8
];


%% load solution from strumpack
nvec=1; % the first rhs is a point source, the second rhs is a gaussian pulse, the third is all-1 inside a square


w=2*pi*f;
wavelength=1/f/smax;
for nth=nvec-1:nvec-1
close all

linefield = {};
xs={};
legends={};
marker={};
Nfd = 0;




%% load solution from BF_VIE_matrix
ppw=wavelength/h1

filename = ['./VIE_F_inc_f_',num2str(f),'_vs_',num2str(vs),'_ivelo_',num2str(ivelo),'_h_',num2str(h1),'_tol_',num2str(tol),'_nth_',num2str(nth),'_matrix.bin'];


fid=fopen(filename);

sizes=fread(fid,3,'int');
h1=fread(fid,1,'double');
tmp=fread(fid,2*sizes(1)*sizes(2)*sizes(3),'double');
M = sizes(1)-1;
N = sizes(2)-1;
K = sizes(3)-1;
F_inc_real=reshape(tmp(1:2:end),M+1,N+1,K+1);
F_inc_imag=reshape(tmp(2:2:end),M+1,N+1,K+1);
F_inc_BF_matrix = F_inc_real + 1i*F_inc_imag;
% F_inc_BF_matrix = permute(F_inc_BF_matrix,[2 3 1]);
F_inc_BF_matrix = F_inc_BF_matrix(1:end-1,1:end-1,1:end-1);  % index end is the right and upper boundary


% if(nth==0)
% F_inc_BF_matrix(floor(M/2)+1,floor(N/2)+1,floor(K/2)+1)=0; % skip source point
% end

step=max(h2s(1),h1)/h1;
F_inc_BF_matrix = F_inc_BF_matrix(1:step:end,1:step:end,1:step:end);

vs1=0;
filename = ['./VIE_F_sca_f_',num2str(f),'_vs_',num2str(vs1),'_ivelo_',num2str(ivelo),'_shape_',num2str(shape),'_h_',num2str(h1),'_tol_',num2str(tol),'_nth_',num2str(nth),'_matrix.bin',twoconst];



fid=fopen(filename);

sizes=fread(fid,3,'int');
h1=fread(fid,1,'double');
tmp=fread(fid,2*sizes(1)*sizes(2)*sizes(3),'double');
M = sizes(1)-1;
N = sizes(2)-1;
K = sizes(3)-1;
F_sca_real=reshape(tmp(1:2:end),M+1,N+1,K+1);
F_sca_imag=reshape(tmp(2:2:end),M+1,N+1,K+1);
F_sca_BF_matrix = F_sca_real + 1i*F_sca_imag;
% F_sca_BF_matrix = permute(F_sca_BF_matrix,[2 3 1]);
F_sca_BF_matrix = F_sca_BF_matrix(1:end-1,1:end-1,1:end-1);  % index end is the right and upper boundary

% 
% if(nth==0)
% F_sca_BF_matrix(floor(M/2)+1,floor(N/2)+1,floor(K/2)+1)=0; % skip source point
% end

step=max(h2s(1),h1)/h1;
F_sca_BF_matrix = F_sca_BF_matrix(1:step:end,1:step:end,1:step:end);
F_tot_BF_matrix = F_sca_BF_matrix + F_inc_BF_matrix;  % VIE uses F_sca_BF instead of -F_sca_BF
VIE_tot_BF_matrix=F_tot_BF_matrix;



x1 = (0:M-1)*h1;
y1 = (0:N-1)*h1;
z1 = (0:K-1)*h1;
x1 = x1(1:step:end);
y1 = y1(1:step:end);
z1 = z1(1:step:end);
[nx,ny,nz] = size(F_inc_BF_matrix);

figure(13)

[X,Y,Z] = meshgrid(y1,z1,x1);
xslice = [floor(ny/2)]*h1;
yslice = [floor(nz/2)]*h1;
zslice = [floor(nx/2)]*h1;

set(gcf,'position',[0 0, 1000 1000]);
set(gcf,'PaperPositionMode','Auto');
slice(X,Y,Z,log(abs(real(-F_inc_BF_matrix))),xslice,yslice,zslice); shading interp; lighting phong; 
hold on;
patch('Faces',faces,'Vertices',vertices,...
'FaceColor',[1, 1, 1],...
'FaceAlpha',  0.3,...
'EdgeColor','k',...
'LineWidth',2);
view(-40,40)
cbar_handle=colorbar('EastOutside');
set(cbar_handle,'position',[0.91 0.42 0.015 0.32])
%     cmin = min([min(min(a_all1(yslice,:,:))), min(min(a_all1(:,xslice,:))), min(min(a_all1(:,:,zslice)))]);
%     cmax = max([max(max(a_all1(yslice,:,:))), max(max(a_all1(:,xslice,:))), max(max(a_all1(:,:,zslice)))]);
%     caxis([cmin,cmax])  
%    caxis([-1,1])
%     set(get(cbar_handle,'ylabel'),'string',FSize,'FontWeight',Fweight);
%     set(get(cbar_handle,'ylabel'),'string','Re(p)','fontsize',FSize,'FontWeight',Fweight);    
gca = get(gcf,'CurrentAxes');
set(get(gca,'title'),'FontSize',FSize,'FontWeight',Fweight,'interpreter','latex');
set(get(gca,'Ylabel'),'FontSize',FSize,'FontWeight',Fweight,'interpreter','latex');
set(get(gca,'Xlabel'),'FontSize',FSize,'FontWeight',Fweight,'interpreter','latex');
set(get(gca,'Zlabel'),'FontSize',FSize,'FontWeight',Fweight,'interpreter','latex');
set(gca,'FontSize',FSize,'FontWeight',Fweight);
set(gca,'Box','on');
set(gca,'Linewidth',1.0);

xlabel('$y$')
ylabel('$x$')
zlabel('$z$')


% colormap hot;
colormap(seismic);
axis equal;
axis ij
grid off



figure(14)

[X,Y,Z] = meshgrid(y1,z1,x1);
xslice = [floor(ny/2)]*h1;
yslice = [floor(nz/2)]*h1;
zslice = [floor(nx/2)]*h1;

set(gcf,'position',[0 0, 1000 1000]);
set(gcf,'PaperPositionMode','Auto');
slice(X,Y,Z,log(abs(real(-F_tot_BF_matrix))),xslice,yslice,zslice); shading interp; lighting phong; 
hold on;
patch('Faces',faces,'Vertices',vertices,...
'FaceColor',[1, 1, 1],...
'FaceAlpha',  0.3,...
'EdgeColor','k',...
'LineWidth',2);
view(-40,40)
cbar_handle=colorbar('EastOutside');
set(cbar_handle,'position',[0.91 0.42 0.015 0.32])
%     cmin = min([min(min(a_all1(yslice,:,:))), min(min(a_all1(:,xslice,:))), min(min(a_all1(:,:,zslice)))]);
%     cmax = max([max(max(a_all1(yslice,:,:))), max(max(a_all1(:,xslice,:))), max(max(a_all1(:,:,zslice)))]);
%     caxis([cmin,cmax])  
%    caxis([-1,1])
%     set(get(cbar_handle,'ylabel'),'string',FSize,'FontWeight',Fweight);
%     set(get(cbar_handle,'ylabel'),'string','Re(p)','fontsize',FSize,'FontWeight',Fweight);    
gca = get(gcf,'CurrentAxes');
set(get(gca,'title'),'FontSize',FSize,'FontWeight',Fweight,'interpreter','latex');
set(get(gca,'Ylabel'),'FontSize',FSize,'FontWeight',Fweight,'interpreter','latex');
set(get(gca,'Xlabel'),'FontSize',FSize,'FontWeight',Fweight,'interpreter','latex');
set(get(gca,'Zlabel'),'FontSize',FSize,'FontWeight',Fweight,'interpreter','latex');
set(gca,'FontSize',FSize,'FontWeight',Fweight);
set(gca,'Box','on');
set(gca,'Linewidth',1.0);

xlabel('$y$')
ylabel('$x$')
zlabel('$z$')


% colormap hot;
colormap(seismic);
axis equal;
axis ij
grid off




idy1=round(yobs/max(h2s(1),h1))+1;
idz1=round(zobs/max(h2s(1),h1))+1;

linefield{Nfd+1} = [-F_tot_BF_matrix(:,idy1(1),idz1(1)),-F_tot_BF_matrix(:,idy1(2),idz1(2))];
legends{Nfd+1} = ['Matrix PPW=',num2str(ppw)];
xs{Nfd+1}=x1;
marker{Nfd+1}='-';
lineColors(Nfd+1,:)=[0,0,0];







%% Load the model for plotting
% filename = ['./VIE_F_inc_f_',num2str(f),'_vs_',num2str(vs),'_ivelo_',num2str(ivelo),'_h_',num2str(h1),'_tol_',num2str(tol),'_nth_',num2str(nth),'_tensor.bin'];
% fid=fopen(filename);
% sizes=fread(fid,3,'int');
% M = sizes(1)-1;
% N = sizes(2)-1;
% K = sizes(3)-1;
% fclose(fid);

x1 = (0:M-1)*h1;
y1 = (0:N-1)*h1;
z1 = (0:K-1)*h1;

if(ivelo==9)
    slowness = zeros(M,N,K);
    slowness(:)=2;
    for ii=1:M
        for jj=1:N
            for kk=1:K
                s0 = 2;
                g1=-0.2; g2=-0.4; g3=-0.35;
                if(x1(ii)<=L && y1(jj)<=L && z1(kk)<=L)
                    slowness(ii,jj,kk) =1/((1.0/s0+g1*(x1(ii)-L/2)+g2*(y1(jj)-L/2)+g3*(z1(kk)-L/2)));
                end
            end
        end
    end
elseif(ivelo==11)
    
    filename = ['slowness_map_shapeoutter',num2str(M+1),'x',num2str(N+1),'x',num2str(K+1),'_off0x0x0_range2_3_nshape20.bin'];
    fid=fopen(filename);
    
    tmp=fread(fid,sizes(1)*sizes(2)*sizes(3),'double');
    slowness=reshape(tmp,M+1,N+1,K+1);
    slowness = slowness(1:end-1,1:end-1,1:end-1);  % index end is the right and upper boundary

end



[nx,ny,nz] = size(slowness);

figure(10)

[X,Y,Z] = meshgrid(y1,z1,x1);
xslice = [floor(ny/2)]*h1;
yslice = [floor(nz/2)]*h1;
zslice = [floor(nx/2)]*h1;

set(gcf,'position',[0 0, 1000 1000]);
set(gcf,'PaperPositionMode','Auto');
slice(X,Y,Z,slowness,xslice,yslice,zslice); shading interp; lighting phong; 
hold on;
patch('Faces',faces,'Vertices',vertices,...
'FaceColor',[1, 1, 1],...
'FaceAlpha',  0.3,...
'EdgeColor','k',...
'LineWidth',2);
view(-40,40)
cbar_handle=colorbar('EastOutside');
set(cbar_handle,'position',[0.91 0.42 0.015 0.32])
gca = get(gcf,'CurrentAxes');
set(get(gca,'title'),'FontSize',FSize,'FontWeight',Fweight,'interpreter','latex');
set(get(gca,'Ylabel'),'FontSize',FSize,'FontWeight',Fweight,'interpreter','latex');
set(get(gca,'Xlabel'),'FontSize',FSize,'FontWeight',Fweight,'interpreter','latex');
set(get(gca,'Zlabel'),'FontSize',FSize,'FontWeight',Fweight,'interpreter','latex');
set(gca,'FontSize',FSize,'FontWeight',Fweight);
set(gca,'Box','on');
set(gca,'Linewidth',1.0);

xlabel('$y$')
ylabel('$x$')
zlabel('$z$')


% colormap hot;
colormap(seismic);
axis equal;
axis ij
grid off















% %% load solution from BF_VIE_tensor
% 
% filename = ['./VIE_F_inc_f_',num2str(f),'_vs_',num2str(vs),'_ivelo_',num2str(ivelo),'_h_',num2str(h1),'_tol_',num2str(tol),'_nth_',num2str(nth),'_tensor.bin'];
% 
% 
% fid=fopen(filename);
% 
% sizes=fread(fid,3,'int');
% h1=fread(fid,1,'double');
% tmp=fread(fid,2*sizes(1)*sizes(2)*sizes(3),'double');
% M = sizes(1)-1;
% N = sizes(2)-1;
% K = sizes(3)-1;
% F_inc_real=reshape(tmp(1:2:end),M+1,N+1,K+1);
% F_inc_imag=reshape(tmp(2:2:end),M+1,N+1,K+1);
% F_inc_BF_tensor = F_inc_real + 1i*F_inc_imag;
% % F_inc_BF_tensor = permute(F_inc_BF_tensor,[2 3 1]);
% F_inc_BF_tensor = F_inc_BF_tensor(1:end-1,1:end-1,1:end-1);  % index end is the right and upper boundary
% 
% 
% % if(nth==0)
% % F_inc_BF_tensor(floor(M/2)+1,floor(N/2)+1,floor(K/2)+1)=0; % skip source point
% % end
% 
% step=max(h2s(1),h1)/h1;
% F_inc_BF_tensor = F_inc_BF_tensor(1:step:end,1:step:end,1:step:end);
% 
% vs1=0;
% filename = ['./VIE_F_sca_f_',num2str(f),'_vs_',num2str(vs1),'_ivelo_',num2str(ivelo),'_shape_',num2str(shape),'_h_',num2str(h1),'_tol_',num2str(tol),'_nth_',num2str(nth),'_tensor.bin',twoconst];
% 
% 
% 
% fid=fopen(filename);
% 
% sizes=fread(fid,3,'int');
% h1=fread(fid,1,'double');
% tmp=fread(fid,2*sizes(1)*sizes(2)*sizes(3),'double');
% M = sizes(1)-1;
% N = sizes(2)-1;
% K = sizes(3)-1;
% F_sca_real=reshape(tmp(1:2:end),M+1,N+1,K+1);
% F_sca_imag=reshape(tmp(2:2:end),M+1,N+1,K+1);
% F_sca_BF_tensor = F_sca_real + 1i*F_sca_imag;
% %F_sca_BF_tensor = permute(F_sca_BF_tensor,[2 3 1]);
% F_sca_BF_tensor = F_sca_BF_tensor(1:end-1,1:end-1,1:end-1);  % index end is the right and upper boundary
% 
% 
% % if(nth==0)
% % F_sca_BF_tensor(floor(M/2)+1,floor(N/2)+1,floor(K/2)+1)=0; % skip source point
% % end
% 
% step=max(h2s(1),h1)/h1;
% F_sca_BF_tensor = F_sca_BF_tensor(1:step:end,1:step:end,1:step:end);
% F_tot_BF_tensor = F_sca_BF_tensor + F_inc_BF_tensor;  % VIE uses F_sca_BF instead of -F_sca_BF
% VIE_tot_BF_tensor=F_tot_BF_tensor;
% 
% 
% 
% x1 = (0:M-1)*h1;
% y1 = (0:N-1)*h1;
% z1 = (0:K-1)*h1;
% x1 = x1(1:step:end);
% y1 = y1(1:step:end);
% z1 = z1(1:step:end);
% [nx,ny,nz] = size(F_inc_BF_tensor);
% 
% figure(17)
% 
% [X,Y,Z] = meshgrid(y1,z1,x1);
% xslice = [floor(ny/2)]*h1;
% yslice = [floor(nz/2)]*h1;
% zslice = [floor(nx/2)]*h1;
% 
% set(gcf,'position',[0 0, 1000 1000]);
% set(gcf,'PaperPositionMode','Auto');
% slice(X,Y,Z,log(abs(real(-F_inc_BF_tensor))),xslice,yslice,zslice); shading interp; lighting phong; 
% hold on;
% patch('Faces',faces,'Vertices',vertices,...
% 'FaceColor',[1, 1, 1],...
% 'FaceAlpha',  0.3,...
% 'EdgeColor','k',...
% 'LineWidth',2);
% view(-40,40)
% cbar_handle=colorbar('EastOutside');
% set(cbar_handle,'position',[0.91 0.42 0.015 0.32])
% %     cmin = min([min(min(a_all1(yslice,:,:))), min(min(a_all1(:,xslice,:))), min(min(a_all1(:,:,zslice)))]);
% %     cmax = max([max(max(a_all1(yslice,:,:))), max(max(a_all1(:,xslice,:))), max(max(a_all1(:,:,zslice)))]);
% %     caxis([cmin,cmax])  
% %    caxis([-1,1])
% %     set(get(cbar_handle,'ylabel'),'string',FSize,'FontWeight',Fweight);
% %     set(get(cbar_handle,'ylabel'),'string','Re(p)','fontsize',FSize,'FontWeight',Fweight);    
% gca = get(gcf,'CurrentAxes');
% set(get(gca,'title'),'FontSize',FSize,'FontWeight',Fweight,'interpreter','latex');
% set(get(gca,'Ylabel'),'FontSize',FSize,'FontWeight',Fweight,'interpreter','latex');
% set(get(gca,'Xlabel'),'FontSize',FSize,'FontWeight',Fweight,'interpreter','latex');
% set(get(gca,'Zlabel'),'FontSize',FSize,'FontWeight',Fweight,'interpreter','latex');
% set(gca,'FontSize',FSize,'FontWeight',Fweight);
% set(gca,'Box','on');
% set(gca,'Linewidth',1.0);
% 
% xlabel('$y$')
% ylabel('$x$')
% zlabel('$z$')
% 
% 
% % colormap hot;
% colormap(seismic);
% axis equal;
% axis ij
% grid off
% 
% 
% 
% figure(18)
% 
% [X,Y,Z] = meshgrid(y1,z1,x1);
% xslice = [floor(ny/2)]*h1;
% yslice = [floor(nz/2)]*h1;
% zslice = [floor(nx/2)]*h1;
% 
% set(gcf,'position',[0 0, 1000 1000]);
% set(gcf,'PaperPositionMode','Auto');
% slice(X,Y,Z,log(abs(real(-F_tot_BF_tensor))),xslice,yslice,zslice); shading interp; lighting phong; 
% hold on;
% patch('Faces',faces,'Vertices',vertices,...
% 'FaceColor',[1, 1, 1],...
% 'FaceAlpha',  0.3,...
% 'EdgeColor','k',...
% 'LineWidth',2);
% view(-40,40)
% cbar_handle=colorbar('EastOutside');
% set(cbar_handle,'position',[0.91 0.42 0.015 0.32])
% %     cmin = min([min(min(a_all1(yslice,:,:))), min(min(a_all1(:,xslice,:))), min(min(a_all1(:,:,zslice)))]);
% %     cmax = max([max(max(a_all1(yslice,:,:))), max(max(a_all1(:,xslice,:))), max(max(a_all1(:,:,zslice)))]);
% %     caxis([cmin,cmax])  
% %    caxis([-1,1])
% %     set(get(cbar_handle,'ylabel'),'string',FSize,'FontWeight',Fweight);
% %     set(get(cbar_handle,'ylabel'),'string','Re(p)','fontsize',FSize,'FontWeight',Fweight);    
% gca = get(gcf,'CurrentAxes');
% set(get(gca,'title'),'FontSize',FSize,'FontWeight',Fweight,'interpreter','latex');
% set(get(gca,'Ylabel'),'FontSize',FSize,'FontWeight',Fweight,'interpreter','latex');
% set(get(gca,'Xlabel'),'FontSize',FSize,'FontWeight',Fweight,'interpreter','latex');
% set(get(gca,'Zlabel'),'FontSize',FSize,'FontWeight',Fweight,'interpreter','latex');
% set(gca,'FontSize',FSize,'FontWeight',Fweight);
% set(gca,'Box','on');
% set(gca,'Linewidth',1.0);
% 
% xlabel('$y$')
% ylabel('$x$')
% zlabel('$z$')
% 
% 
% % colormap hot;
% colormap(seismic);
% axis equal;
% axis ij
% grid off
% 
% 
% 
% 
% idy1=round(yobs/max(h2s(1),h1))+1;
% idz1=round(zobs/max(h2s(1),h1))+1;
% 
% linefield{Nfd+2} = [-F_tot_BF_tensor(:,idy1(1),idz1(1)),-F_tot_BF_tensor(:,idy1(2),idz1(2))];
% legends{Nfd+2} = ['Tensor PPW=',num2str(ppw)];
% xs{Nfd+2}=x1;
% marker{Nfd+2}='-';
% lineColors(Nfd+2,:)=[0,0,0];




% axisticksize = 40;
% origin = [200,60];
% markersize = 6;
% LineWidth = 2;
% 
% 
% figure(15)
% legends={};
% for nn=1:size(linefield{2},2)
%     
%     plot(xs{2},(log10(abs(real(linefield{2}(:,nn))))),marker{2},'LineWidth',LineWidth,'Color',[lineColors(nn,:)],'MarkerFaceColor',[lineColors(nn,:)],'MarkerSize',markersize)
%     hold on
%     xlim([min(xs{2}),max(xs{2})])
%     
%     str=['y: ',num2str(yobs(nn)), ' z: ',num2str(zobs(nn))];
%     legends{nn}=str;
% end
% 
% xlabel('$x$')
% ylabel('log$|\mathrm{real}(u)|$')
% 
% legend(legends,'NumColumns',1,'Location','SouthEast')
% legend('boxoff')
% gca = get(gcf,'CurrentAxes');
% set( gca, 'FontName','Times New Roman','fontsize',axisticksize);
% 
% set(get(gca,'Ylabel'),'FontSize',FSize,'FontWeight',Fweight,'interpreter','latex');
% set(get(gca,'Xlabel'),'FontSize',FSize,'FontWeight',Fweight,'interpreter','latex');
% 
% 
% set(gcf,'Position',[origin,1000,700]);


% cmin = min(min(min(min(((real(a_all0)))))),min(min(min(((real(-F_inc_BF_matrix)))))));
% cmax = max(max(max(max(((real(a_all0)))))),max(max(max(((real(-F_inc_BF_matrix)))))));


% cmin = min(min(min(((real(-F_inc_BF_matrix))))));
% cmax = max(max(max(((real(-F_inc_BF_matrix))))));


figure(13)
saveas(gcf,['Matrix_wavefield_3D_inc_f_',num2str(f),'_ivelo_',num2str(ivelo),'_h_',num2str(h2s(1)),'_nth_',num2str(nth),'_ppw',num2str(wavelength/h1),'BF.png'])    
figure(14)
saveas(gcf,['Matrix_wavefield_3D_tot_f_',num2str(f),'_ivelo_',num2str(ivelo),'_h_',num2str(h2s(1)),'_nth_',num2str(nth),'_ppw',num2str(wavelength/h1),'BF.png'])    
% figure(17)
% saveas(gcf,['Tensor_wavefield_3D_inc_f_',num2str(f),'_ivelo_',num2str(ivelo),'_h_',num2str(h2s(1)),'_nth_',num2str(nth),'_ppw',num2str(wavelength/h1),'BF.png'])    
% figure(18)
% saveas(gcf,['Tensor_wavefield_3D_tot_f_',num2str(f),'_ivelo_',num2str(ivelo),'_h_',num2str(h2s(1)),'_nth_',num2str(nth),'_ppw',num2str(wavelength/h1),'BF.png'])    
% 
% 
% 
% figure(15)
% saveas(gcf,['Linefield_inc_3D_f_',num2str(f),'_ivelo_',num2str(ivelo),'_nth_',num2str(nth),'_line_1_ppw',num2str(wavelength/h1),'BF_Tensoronly.png'])    


figure(10)
saveas(gcf,['Slowness_f_',num2str(f),'_ivelo_',num2str(ivelo),'_h_',num2str(h2s(1)),'.png'])    


end


