function K = covSEiso_noise(hyp, x, z, i)

% Squared Exponential covariance function with isotropic distance measure. The
% covariance function is parameterized as:
%
% k(x^p,x^q) = sf^2 * exp(-(x^p - x^q)'*inv(P)*(x^p - x^q)/2) 
%
% where the P matrix is ell^2 times the unit matrix and sf^2 is the signal
% variance. The hyperparameters are:
%
% hyp = [ log(ell)
%         log(sf)  ]
%
% For more help on design of covariance functions, try "help covFunctions".
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2010-09-10.
%
% See also COVFUNCTIONS.M.

if nargin<2, K = '3'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = isempty(z); dg = strcmp(z,'diag');                       % determine mode

ell = exp(hyp(1));                                 % characteristic length scale
sf2 = exp(2*hyp(2));                                           % signal variance
% noise = hyp(3);

% precompute squared distances
if dg                                                               % vector kxx
  K = zeros(size(x,1),1);
else
  if xeqz                                                 % symmetric matrix Kxx
    K = sq_dist(x'/ell) + add_noise(hyp, x', x');
  else                                                   % cross covariances Kxz
    K = sq_dist(x'/ell,z'/ell)+ add_noise(hyp, x', z');
  end
end

if nargin<4                                                        % covariances
  K = sf2*exp(-K/2);
else                                                               % derivatives
  if i==1
    K = sf2*exp(-K/2).*K;
  elseif i==2
    K = 2*sf2*exp(-K/2);
  else
    error('Unknown hyperparameter')
  end
end

end

function f = add_noise(hyp, xa, xb)
f = zeros(size(xa,2), size(xb,2));
ell = exp(hyp(1));                                 % characteristic length scale
sf2 = exp(2*hyp(2));                                           % signal variance
noise = hyp(3);
% global in_noise
% noise=in_noise;
for i=1:size(xa,2),
    XA = repmat(xa(:,i),1,size(xb,2));
    d = XA - xb;
    s = sum(d.*d, 1);
    
    %%calculating double derivative of covariance function
    w = 1/ell^2;
    temp = w * sf2 * exp(-s/(2*ell^2));
    doubleDerCov = temp .*(w*s - 3);

    %%calculating quad derivative of covariance function
    
    quadDerCovariance = quadDerCov(hyp, XA, xb, s, w);
    
    % gp.s and gp.l(hyperparameters) govern properties of sample function
    % gp.s controls the typical amplitude
    % gp.l controls the lengthscale of the variation
   f(i,:) = noise*doubleDerCov + noise^2/4*quadDerCovariance;
end

end

function f = quadDerCov(hyp, xa, xb, s, w)
    ell = exp(hyp(1));                                 % characteristic length scale
    sf2 = exp(2*hyp(2));                               % signal variance

    temp1 = w^2 * sf2 * exp(-s/(2*ell^2));
    temp2 = sf2 * w^2 * (9 + w*( w * xa(1,:).^4 + w * xa(2,:).^4 -4*w * xa(1,:).^3.*xb(1,:) - 6*(xa(3,:).^2 + xb(1,:).^2) - 4*xa(1,:).*xb(1,:) ...
        .*(-3 + w*xb(1,:).^2) + 6*xa(1,:).^2 .* (-1 + w*xb(1,:).^2) - 4*w*xa(2,:).^3.*xb(2,:) - 6*xb(2,:).^2 - ...
         4*xa(2,:).*xb(2,:).*(-3+w*xb(2,:).^2) + 6*(-1+w*xa(3,:).^2).*xb(3,:).^2  - 4*w*xa(3,:).*xb(3,:).^3 + w*xb(3,:).^4 ));
    f = temp1 .* temp2;
end
