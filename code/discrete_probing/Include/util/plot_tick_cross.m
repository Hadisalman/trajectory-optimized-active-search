function [h] = plot_tick_cross(offset, scale,level,alpha,beta,reward,direction)

h = [];

if reward == 1
  clr = 'r';
  x_shape = [-1  -.2  1];
  y_shape = [.2   -1   1];
else
  clr = 'g';%[.7 .7 .7];
  x_shape = [-1 -1; 1  1];
  y_shape = [-1  1; 1 -1];
end

x_target = offset(1) + .5*scale*x_shape;
y_target = offset(2) + .5*scale*y_shape; 

h = [];

h = [h, plot( x_target, y_target,'color',clr,'linewidth',1)];

