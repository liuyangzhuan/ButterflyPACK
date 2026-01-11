close all;
clear all;
clc;
format shortE
FSize=50;
Fweight='normal';
load './seismic.map'
twoconst=[''];


lineColors = line_colors(6);



% 
h2s = [0.00125];
% yobs = [0.6 0.52-max(h2s)*10];
% zobs = [0.6 0.52-max(h2s)*10];
% yobs = [0.6 max(h2s)*10];
% zobs = [0.6 max(h2s)*10];

f=40;
vs=1;
ivelo=1;
smax=2;
tol=1e-2;  % tolerence in weno_bf
h1=h2s(1);  %0.015625; %0.015625; % grid size in weno_bf
h3=h2s(1);  % grid size in 9-point fdfd
L=0.5;
center=[0.25,0.25,0.25];




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
for nth=0:nvec-1
close all

linefield = {};
xs={};
legends={};
marker={};
Nfd = 0;




%% load solution from BF_VIE_matrix
ppw=wavelength/h1;


%% Load the model for plotting
filename = ['../build/F_inc_f_',num2str(f),'_vs_',num2str(vs),'_ivelo_',num2str(ivelo),'_h_',num2str(h1),'_tol_',num2str(tol),'_nth_',num2str(nth),'_tensor.bin'];
fid=fopen(filename);
sizes=fread(fid,3,'int');
M = sizes(1)-1;
N = sizes(2)-1;
K = sizes(3)-1;
fclose(fid);

x1 = (0:M-1)*h1;
y1 = (0:N-1)*h1;
z1 = (0:K-1)*h1;


fid=fopen(filename);

sizes=fread(fid,3,'int');
h1=fread(fid,1,'double');
tmp=fread(fid,2*sizes(1)*sizes(2)*sizes(3),'double');
M = sizes(1)-1;
N = sizes(2)-1;
K = sizes(3)-1;
F_inc_real=reshape(tmp(1:2:end),M+1,N+1,K+1);
F_inc_imag=reshape(tmp(2:2:end),M+1,N+1,K+1);
F_inc_BF_tensor = F_inc_real + 1i*F_inc_imag;
% F_inc_BF_tensor = permute(F_inc_BF_tensor,[2 3 1]);
F_inc_BF_tensor = F_inc_BF_tensor(1:end-1,1:end-1,1:end-1);  % index end is the right and upper boundary


% if(nth==0)
% F_inc_BF_tensor(floor(M/2)+1,floor(N/2)+1,floor(K/2)+1)=0; % skip source point
% end

step=max(h2s(1),h1)/h1;
F_inc_BF_tensor = F_inc_BF_tensor(1:step:end,1:step:end,1:step:end);




x1 = (0:M-1)*h1;
y1 = (0:N-1)*h1;
z1 = (0:K-1)*h1;
x1 = x1(1:step:end);
y1 = y1(1:step:end);
z1 = z1(1:step:end);
[nx,ny,nz] = size(F_inc_BF_tensor);

figure(1)

[X,Y,Z] = meshgrid(y1,z1,x1);
xslice = [floor(ny/2)]*h1;
yslice = [floor(nz/2)]*h1;
zslice = [floor(nx/2)]*h1;

set(gcf,'position',[0 0, 1000 1000]);
set(gcf,'PaperPositionMode','Auto');
slice(X,Y,Z,log(abs(real(-F_inc_BF_tensor))),xslice,yslice,zslice); shading interp; lighting phong; 
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


figure(1)
saveas(gcf,['Tensor_wavefield_3D_inc_f_',num2str(f),'_ivelo_',num2str(ivelo),'_h_',num2str(h2s(1)),'_nth_',num2str(nth),'_BF.png'])    


end


