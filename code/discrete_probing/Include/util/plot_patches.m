function [] = plot_patches(regions, vals, varargin)

hold on;

if nargin<2
  
  colors = {[166 206 227], [31 120 180], [178 223 138], [51 160 44], [251 154 153], ...
    [227 26 28], [253 191 111], [255 127 0], [202 178 214], [106 61 154], ...
    [255 255 153]};
  
  for i=1:size(regions,1)
    reg_x = regions(i,1:2);
    reg_y = regions(i,3:4);
    patch([reg_x(1) reg_x(2) reg_x(2) reg_x(1)], ...
      [reg_y(1) reg_y(1) reg_y(2) reg_y(2)], colors{i}/255);
  end

elseif strcmpi(vals, 'rule')
  
  for i=1:size(regions,1)
    reg_x = regions(i,1:2);
    reg_y = regions(i,3:4);
    plot([reg_x(1) reg_x(2) reg_x(2) reg_x(1), reg_x(1)], ...
      [reg_y(1) reg_y(1) reg_y(2) reg_y(2), reg_y(1)]  , 'b'    );
  end

else
  for i=1:size(regions,1)
    reg_x = regions(i,1:2);
    reg_y = regions(i,3:4);
    patch([reg_x(1) reg_x(2) reg_x(2) reg_x(1)], ...
      [reg_y(1) reg_y(1) reg_y(2) reg_y(2)], colorlookup(vals(i), varargin{:}));
  end

end
