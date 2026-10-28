function plot_tensor_VIE_slowness(slowness_file, n, h, outprefix, fracs)
% Slices of the slowness model of an ivelo 11 tensor VIE run (the file of EXAMPLE/slowness_model3d_generator.py), in
% the coordinates the driver uses, so that they line up with plot_tensor_VIE_slices.m's field planes:
%  - the three centre planes in 3D (same view and box as the field figures);
%  - for each of x, y and z, the planes at the fractions `fracs` of the box (default 1/8, 3/8, 1/2, 5/8, 7/8).
% The driver reads value (i,j,k) (0-based, x fastest) at offset i + I*j + I*J*k: MATLAB's column-major [I J K] array.
% (numpy wrote it with its last axis fastest, so the generator's own plot shows these axes in the reverse order.)
% n = [I J K] grid points, h = grid spacing.  The file is memory-mapped: only the slices are read.
% Example: plot_tensor_VIE_slowness('slowness_map_..._nshape25.bin', [320 320 320], 0.00125, 'n320')
if nargin < 5, fracs = [1 3 4 5 7] / 8; end
m = memmapfile(slowness_file, 'Format', {'double', double(n(:))', 's'});
I = n(1); J = n(2); K = n(3);
x = (0:I-1)*h;  y = (0:J-1)*h;  z = (0:K-1)*h;
ci = floor(I/2); cj = floor(J/2); ck = floor(K/2);   % the field files' centre indices (1-based: n/2)
Sx = squeeze(m.Data.s(ci,:,:));  Sy = squeeze(m.Data.s(:,cj,:));  Sz = m.Data.s(:,:,ck);
% the planes along each axis, read first so that one colour range covers every figure
nf = numel(fracs);
P = cell(3, nf);  pos = zeros(3, nf);
for a = 1:3
  for f = 1:nf
    idx = max(1, min(n(a), round(fracs(f)*n(a))));
    pos(a,f) = (idx-1)*h;
    switch a
      case 1, P{a,f} = squeeze(m.Data.s(idx,:,:));
      case 2, P{a,f} = squeeze(m.Data.s(:,idx,:));
      case 3, P{a,f} = m.Data.s(:,:,idx);
    end
  end
end
vals = [Sx(:); Sy(:); Sz(:); cell2mat(cellfun(@(q) q(:), P(:), 'UniformOutput', false))];
smin = min(vals);  smax = max(vals);
FSize = 20;

% 3D view of the centre planes
figure('visible','off','position',[0 0 1000 1000]);
[Yg,Zg] = ndgrid(y,z);  surf(x(ci)*ones(size(Yg)), Yg, Zg, Sx, 'EdgeColor','none'); hold on;
[Xg,Zg] = ndgrid(x,z);  surf(Xg, y(cj)*ones(size(Xg)), Zg, Sy, 'EdgeColor','none');
[Xg,Yg] = ndgrid(x,y);  surf(Xg, Yg, z(ck)*ones(size(Xg)), Sz, 'EdgeColor','none');
L = [x(end) y(end) z(end)];
v = [0 0 0; L(1) 0 0; L(1) L(2) 0; 0 L(2) 0; 0 0 L(3); L(1) 0 L(3); L(1) L(2) L(3); 0 L(2) L(3)];
patch('Faces',[1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8],'Vertices',v,'FaceColor',[1 1 1], ...
      'FaceAlpha',0.15,'EdgeColor','k','LineWidth',2);
view(-40,40); colormap(parula); caxis([smin smax]); cb = colorbar('EastOutside'); ylabel(cb,'slowness');
axis equal; grid off; box on;
xlabel('$x$','interpreter','latex'); ylabel('$y$','interpreter','latex'); zlabel('$z$','interpreter','latex');
title({'slowness, centre planes', sprintf('%.3g .. %.3g', smin, smax)},'FontSize',14);
set(gca,'FontSize',FSize);
saveas(gcf, sprintf('%s_slowness3d.png', outprefix));
close(gcf);

% planes at several positions along each axis
figure('visible','off','position',[0 0 380*nf 1150]);
names = {'x','y','z'};
for a = 1:3
  switch a
    case 1, u = y; w = z; lu = 'y'; lw = 'z';
    case 2, u = x; w = z; lu = 'x'; lw = 'z';
    case 3, u = x; w = y; lu = 'x'; lw = 'y';
  end
  for f = 1:nf
    subplot(3,nf,(a-1)*nf+f);
    imagesc(u, w, P{a,f}.'); axis xy equal tight; colormap(gca,parula); caxis([smin smax]);
    colorbar;
    title(sprintf('%s = %.3f', names{a}, pos(a,f)),'FontSize',12); xlabel(lu); ylabel(lw);
  end
end
saveas(gcf, sprintf('%s_slowness_planes.png', outprefix));
close(gcf);
end
