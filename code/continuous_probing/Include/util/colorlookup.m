function [rgb] = colorlookup(val, clim, cmap)

if nargin<2, clim = get(gca,'clim'); end

if nargin<3, cmap = colormap; end

ind = round( ...
  (val - clim(1)) / (clim(2)-clim(1)) * size(cmap,1) ...
  );

ind = max(ind, 1);
ind = min(ind, size(cmap,1));

rgb = cmap(ind, :);

end
