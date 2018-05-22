function [] = plot_demo(x_shape, x_view, f_view, fs2_view, x_collected, regions, level, alpha, beta2, found, varargin)

p = inputParser;
p.addParamValue('clims',            [-1, 3]);
p.addParamValue('plotTailGaussian', true);
p.addParamValue('fs2_lims',         [0, 3]);

p.parse(varargin{:});
clims            = p.Results.clims;
plotTailGaussian = p.Results.plotTailGaussian;
fs2_lims         = p.Results.fs2_lims;

granularity = regions(1,2) - regions(1,1);
siz         = [min(regions(:,1)), max(regions(:,2)), min(regions(:,3)), max(regions(:,4))];

% ------------------- begin plotting -----------------------

image([siz(1), siz(2)], [siz(3),siz(4)], ...
  reshape(colorlookup(sqrt(max(fs2_view,0)), fs2_lims, 1-gray), ...
  [x_shape,3] ));
set(gca,'yDir','normal');
hold on;

caxis(clims); colormap(brighten(jet, .7));
plot_patches(regions, 'rule'); hold on;

% [c,h] = contour( reshape(x_view(:,1), x_shape),  reshape(x_view(:,2), x_shape),  reshape(f_view, x_shape), linspace(clims(1), clims(2), 9));
[c,h] = contour( reshape(x_view(:,1), x_shape),  reshape(x_view(:,2), x_shape),  reshape(f_view, x_shape), [0 0.75 1 1.25]);

set(h,'linewidth',3);

clabel(c,h,'fontsize',20); 

for i=1:length(regions)
  if plotTailGaussian
    plot_tail_gaussian([
      .5*regions(i,1)+.5*regions(i,2)
      .5*regions(i,3)+.5*regions(i,4)
      ], .2*granularity, level, alpha(i), sqrt(beta2(i)), found(i), 1);
  else
    plot_tick_cross([
      .5*regions(i,1)+.5*regions(i,2)
      .5*regions(i,3)+.5*regions(i,4)
      ], .2*granularity, level, alpha(i), sqrt(beta2(i)), found(i), 1);
  end
end

if ~isempty(x_collected)
  scatter(x_collected(:,1), x_collected(:,2), 50, 'k', 's', 'filled');
end

axis([siz(1) siz(2) siz(3) siz(4)]);
drawnow;
