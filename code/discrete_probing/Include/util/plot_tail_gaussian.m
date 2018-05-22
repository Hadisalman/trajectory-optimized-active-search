function [h] = plot_tail_gaussian(offset, scale,level,alpha,beta,reward,direction)

if nargin < 7, direction = 1; end

if reward
  clr = 'r';
else
  clr = 'g';%[.7 .7 .7];
end

% default
x = level + (-3:.5:3)*beta;
y = normpdf(x,alpha,beta);

% trim
% if nargin<7 || isempty(granularity)
%   granularity = 1;
% end
% y = min(y, granularity/scale);

ell = zeros(size(x));

% center level
x_level = x - level;

x_target = offset(1) + scale*x_level;
y_target = offset(2) + scale*y;
ell_target = offset(2) + scale*ell;


h = [];

h(end+1) = plot( x_target, y_target,'color',clr,'linewidth',3);
h(end+1) = plot( x_target, ell_target,'color',clr,'linewidth',3);

h(end+1) = plot( x_target(x==level) * [1 1], ...
                 [y_target(x==level), ell_target(x==level) - .2*scale], ...
                 'color', clr,'linewidth',3);
% h(end+1) = quiver(x_target(x==level), ell_target(x==level) - .2*scale, ...
%   sign(direction)*scale, 0, 'color', clr,'linewidth',3);
h(end+1) = plot( x_target(x==level) + [0, sign(direction)*scale], ...
                 [ell_target(x==level) - .2*scale] + [0 0], ...
                 'color', clr,'linewidth',3);

% if sign(direction)==1
%   h(end+1) = area(x_target(x>=level), y_target(x>=level), ell_target(x==level), ...
%     'facecolor',clr,'edgecolor',clr);
% else
%   h(end+1) = area(x_target(x<=level), y_target(x<=level), ell_target(x==level), ...
%     'facecolor',clr,'edgecolor',clr);
% end

