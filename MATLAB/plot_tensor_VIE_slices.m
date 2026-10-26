function plot_tensor_VIE_slices(sca_files, inc_files, labels, outprefix)
% Plots of the tensor VIE solution from the three centre slices the driver writes (*_tensor_slices.bin, for
% grids too large for the full tensor file), for one or several runs of the same problem:
%  - per run, the three centre planes of log|Re(-F_tot)| in 3D (as figures 17/18 of plot_tensor_VIE_field.m,
%    whose slice() needs the whole volume: ~15 GB a field at 640^3; the planes are drawn directly instead);
%  - Re(F_sca) on the z centre plane of every run side by side and, with two or more runs, each run's
%    |F_sca - F_sca(last run)| on that plane with its relative difference (the last run is the reference).
% sca_files, inc_files: cell arrays of scattered/incident slice files (one incident file is used for all runs
% when only one is given); labels: cell array of titles; outprefix: prefix of the PNG files written.
% Example (640^3, the tensor compressed to 1e-2 and to 1e-3):
%   plot_tensor_VIE_slices({'sca_tc2.bin','sca_tc3.bin'}, {'inc_tc2.bin','inc_tc3.bin'}, ...
%                          {'tol\_comp 1e-2','tol\_comp 1e-3'}, 'n640')
% Run without a display: matlab -batch "addpath('<repo>/MATLAB'); plot_tensor_VIE_slices(...)"
here = fileparts(mfilename('fullpath'));
seismic = load(fullfile(here, 'seismic.map'));
FSize = 20;
nr = numel(sca_files);
S = cell(nr,1);  T = cell(nr,1);
for r = 1:nr
  S{r} = read_slices(sca_files{r});
  I = read_slices(inc_files{min(r, numel(inc_files))});
  T{r} = S{r};
  T{r}.x = S{r}.x + I.x;  T{r}.y = S{r}.y + I.y;  T{r}.z = S{r}.z + I.z;   % VIE uses F_sca instead of -F_sca
end

% the scatterer box: the whole grid
h = S{1}.h;  L = (S{1}.n - 1) * h;
v = [0 0 0; L(1) 0 0; L(1) L(2) 0; 0 L(2) 0; 0 0 L(3); L(1) 0 L(3); L(1) L(2) L(3); 0 L(2) L(3)];
f = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];

for r = 1:nr
  figure('visible','off','position',[0 0 1000 1000]);
  plot_planes(T{r}, @(F) log(abs(real(-F))));
  caxis(log_range(T{r}));
  hold on;
  patch('Faces',f,'Vertices',v,'FaceColor',[1 1 1],'FaceAlpha',0.3,'EdgeColor','k','LineWidth',2);
  view(-40,40); colormap(seismic); colorbar('EastOutside'); axis equal; grid off; box on;
  xlabel('$x$','interpreter','latex'); ylabel('$y$','interpreter','latex'); zlabel('$z$','interpreter','latex');
  title({labels{r}, 'log|Re(-F_{tot})|'},'FontSize',14);
  set(gca,'FontSize',FSize);
  saveas(gcf, sprintf('%s_tot3d_%d.png', outprefix, r));
  close(gcf);
end

rows = 1 + (nr > 1);
figure('visible','off','position',[0 0 560*nr 520*rows]);
x = (0:S{1}.n(1)-1)*h;  y = (0:S{1}.n(2)-1)*h;  zc = (S{1}.idx(3)-1)*h;
ref = S{nr}.z;
clim = max(abs(real(ref(:))));
for r = 1:nr
  subplot(rows,nr,r);
  imagesc(x,y,real(S{r}.z).'); axis xy equal tight; colormap(gca,seismic); caxis([-clim clim]); colorbar;
  title({labels{r}, sprintf('Re(F_{sca}), z = %.3g', zc)},'FontSize',12); xlabel('x'); ylabel('y');
  if nr > 1 && r == nr
    subplot(rows,nr,nr+r); axis off;
    text(0.5, 0.5, 'reference for the differences', 'HorizontalAlignment','center', 'FontSize',12);
  elseif nr > 1
    subplot(rows,nr,nr+r);
    imagesc(x,y,abs(S{r}.z - ref).'); axis xy equal tight; colormap(gca,parula); caxis([0 clim]); colorbar;
    title({'|F_{sca} - F_{sca}(last)|', sprintf('relative %.3f on the plane', norm(S{r}.z(:)-ref(:))/norm(ref(:)))},'FontSize',12);
    xlabel('x'); ylabel('y');
  end
end
saveas(gcf, sprintf('%s_zslice_compare.png', outprefix));
close(gcf);
end

function plot_planes(S, fn)
% the three centre planes of a slice file as surfaces at their physical positions, every second point
% (software OpenGL renders 3 x 640^2 faces a figure slowly)
x = (0:S.n(1)-1)*S.h;  y = (0:S.n(2)-1)*S.h;  z = (0:S.n(3)-1)*S.h;
k1 = 1:2:numel(x);  k2 = 1:2:numel(y);  k3 = 1:2:numel(z);
[Yg,Zg] = ndgrid(y(k2),z(k3));  surf(x(S.idx(1))*ones(size(Yg)), Yg, Zg, fn(S.x(k2,k3)), 'EdgeColor','none'); hold on;
[Xg,Zg] = ndgrid(x(k1),z(k3));  surf(Xg, y(S.idx(2))*ones(size(Xg)), Zg, fn(S.y(k1,k3)), 'EdgeColor','none');
[Xg,Yg] = ndgrid(x(k1),y(k2));  surf(Xg, Yg, z(S.idx(3))*ones(size(Xg)), fn(S.z(k1,k2)), 'EdgeColor','none');
shading interp; lighting phong;
end

function c = log_range(S)
% colour range of log|Re(-F)| over the three planes: 2nd percentile to the maximum (zero crossings give -Inf..-16)
v = log(abs(real(-[S.x(:); S.y(:); S.z(:)])));
v = sort(v(isfinite(v)));
c = [v(max(1,round(0.02*numel(v)))) v(end)];
end

function S = read_slices(filename)
% read_vie_tensor_or_slices of plot_tensor_VIE_field.m, keeping the three slices as 2D arrays
fid = fopen(filename);
if fid < 0, error(['Cannot open ',filename]); end
magic = fread(fid,1,'int'); version = fread(fid,1,'int');
if magic ~= 20260605 || version ~= 1, fclose(fid); error(['Unsupported tensor slice file: ',filename]); end
sizes = fread(fid,3,'int'); idx = fread(fid,3,'int'); h1 = fread(fid,1,'double');
nx = sizes(1); ny = sizes(2); nz = sizes(3);
tmp = fread(fid, 2*(ny*nz + nx*nz + nx*ny), 'double'); fclose(fid);
vals = tmp(1:2:end) + 1i*tmp(2:2:end);
p0 = 0;   p1 = ny*nz;       S.x = reshape(vals(p0+1:p1),ny,nz);
p0 = p1;  p1 = p0 + nx*nz;  S.y = reshape(vals(p0+1:p1),nx,nz);
p0 = p1;  p1 = p0 + nx*ny;  S.z = reshape(vals(p0+1:p1),nx,ny);
S.idx = idx;  S.n = sizes;  S.h = h1;
end
