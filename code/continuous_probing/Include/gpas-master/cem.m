function [x, c, mu, C] = cem(fun, x0, opts, varargin)
% The cross-entropy method
% @param fun function to  be minimized
% @param x0 initial guess
% options:
% @param opts.N: number of samples
% @param opts.rho: quantile (e.g. 0.1)
% @param opts.C: initial covariance
% @param opts.iter: total iterations
% @param opts.v: update coefficients
% @param varagin any other arguments that will be passed to fun
%
% @return x best sample
% @return c best cost
% @return mu distribution mean
% @return C distribution covariance
%
% Author: Marin Kobilarov, marin@jhu.edu

d = length(x0);

if ~isfield(opts, 'N')
  opts.N = d*20;
end

if ~isfield(opts, 'rho')
  opts.rho = 0.1;
end

if ~isfield(opts, 'C')
  opts.C = eye(d);
end

if ~isfield(opts, 'iter')
  fun = 20;
end

if ~isfield(opts, 'v')
  opts.v = 1;
end

if ~isfield(opts, 'sigma')
  opts.sigma = 0;
end

if opts.sigma 
  opts.N = 2*d+1;
end

if ~isfield(opts, 'tilt')
  opts.tilt = 0;
end

if ~isfield(opts, 'lb')
  opts.lb = [];
end

if ~isfield(opts, 'ub')
  opts.ub = [];
end


if opts.tilt 
end


N = opts.N;
nf = round(opts.rho*N);
C = opts.C;
v = opts.v;

cs = zeros(N, 1);
xs = zeros(d, N);

x = x0;
c = inf;

mu = x0;

a = 0.001;
k = 0;
b = 2;
l = a*a*(d+k)-d;
Ws = [l/(d+l), repmat(1/(2*(d+l)), 1, 2*d)];
Wc = [l/(d+l) + (1-a*a+b), repmat(1/(2*(d+l)), 1, 2*d)];


for j=1:opts.iter
  
  if opts.sigma    
    % this is an experimental version of the CE method using sigma-points
    
    A = sqrt(d+l)*chol(C)';
    xs = [mu, repmat(mu, 1, d) + A, repmat(mu, 1, d) - A];

    xm = zeros(d,1);
    for i=1:size(xs,2),
      fi = fun(xs(:,i), varargin{:});
      if (length(fi) > 1)
        cs(i) = sum(fi.*fi);
      else
        cs(i) = fi;
      end      
      cs(i) = exp(-cs(i));
      
      xm = xm + Ws(i)*cs(i)*xs(:,i);      
    end
    
    Pm = zeros(d,d);
    for i=1:size(xs,2),      
      dx = xs(:,i) - xm;
      Pm = Pm + Wc(i)*cs(i)*dx*dx';
    end
    
    csn = sum(cs);
    mu = mu/csn;
    C = Pm/csn;
    
    x = mu;
    c = cs(1);
    
  else
    % this is the standard CE method using random sampling
    
    
    if (~isempty(opts.lb))
      n = length(mu);
      A=[-eye(n); 
         eye(n)];
      B=[-opts.lb; 
         opts.ub];
      xs = rmvnrnd(mu, C, N, A, B)';
    else
      xs = mvnrnd(mu, C, N)';
    end
    
    
    for i=1:N,
      fi = fun(xs(:,i), varargin{:});
      if (length(fi) > 1)
        cs(i) = sum(fi.*fi)/2;
      else
        cs(i) = fi;
      end
    end    
    
    if ~opts.tilt 
      [cs,is] = sort(cs, 1, 'ascend');
      xes = xs(:, is(1:nf));
      
      mu = (1 - v).*mu + v.*mean(xes')';
      C = (1 - v).*C + v.*cov(xes');% + diag(opts.rf.*rand(opts.n,1));    
      
      if (cs(1) < c)
        x = xes(:,1);
        c = cs(1);
      end
    
    else

      if (j==1)
        S.ps0 = mvnpdf(xs', mu', C);
      end            
      
      [cmin, imin] = min(cs);

            
%      b = max(1/cmin, .001);
%      b = max(1/(max(cs)-min(cs)), .001);

      %good one:
      b = 1/mean(cs);
%      b = max(1/min(cs), .001);

      if 0
      b = b*(entropy(mu, C));

      S.ps = mvnpdf(xs',mu', C);
      
      S.Jh = mean(cs);
      S.Js = cs;
      bmin = 0;
      bmax = 1;

      S.xs = xs;
      S.v = v;
      S.mu = mu;
      S.C = C;

      bs = bmin:.001:bmax;
      gs = zeros(size(bs));
      for l=1:length(bs)
        gs(l) = minb3(bs(l), S);
      end

      plot(bs, gs,'g');
      drawnow
      gs
      
      [gm,bi]=min(gs);
      b=bs(bi)
      
      [b,FVAL,EXITFLAG,OUTPUT] = fminbnd(@(b) minb3(b, S), bmin, bmax)      

      keyboard


      global ws
      end
      
      %      kl = sum(-log(ws))/N
      
      %      b = b*kl;
      

      %     b = 1;
      ws = exp(-b*cs);
      ws = ws/sum(ws);
      
      mu = (1 - v).*mu + v.*(xs*ws);
      C = (1 - v).*C + v.*weightedcov(xs', ws);
      
      if (cmin < c)
        x = xs(:,imin);
        c = cmin;
      end
      
    end
    
  end
  
end

function f = minb(b, S)

b
N = length(S.Js);
ws = exp(-b*S.Js);

eta = sum(ws)/N;

delta = .1;
g = sqrt(log(1/delta)/(2*N));

ws = ws/sum(ws);

mu = (1 - S.v).*S.mu + S.v.*(S.xs*ws);

C = (1 - S.v).*S.C + S.v.*weightedcov(S.xs', ws);

C
vs = 1/eta*ws.*(-log(eta)*ones(N,1) - b*S.Js + log(S.ps) ...
                - log(mvnpdf(S.xs', mu', C)));

R = max(vs)-min(vs)

f = sum(vs)/N + R*g;

function f = minb2(b, S)

N = length(S.Js);
ws = exp(-b*S.Js);

eta = sum(ws)/N;

delta = .1;
g = sqrt(log(1/delta)/(2*N));

ws = ws/sum(ws);

mu = (1 - S.v).*S.mu + S.v.*(S.xs*ws);

C = (1 - S.v).*S.C + S.v.*weightedcov(S.xs', ws);

mu 
C
Ws = S.ps0./mvnpdf(S.xs', mu', C);
Ws

vs = exp(-2*b*S.Js).*Ws;

R = max(vs);

f = sum(vs)/N + R*g;



function f = minb3(b, S)

N = length(S.Js);
Jmin = min(S.Js)
Jmax = max(S.Js)

ws = exp(-b*S.Js/Jmax);

delta = .5;
g = sqrt(log(1/delta)/(2*N))

mean(ws) - exp(-b*Jmin/Jmax)*g
b*(mean(S.Js)/Jmax - g)

f = log(mean(ws) - exp(-b*Jmin/Jmax)*g) + b*(mean(S.Js)/Jmax - g);
f
f = -f;


function C = weightedcov(Y, w)
%   Weighted Covariance Matrix
%
%   WEIGHTEDCOV returns a symmetric matrix C of weighted covariances
%   calculated from an input T-by-N matrix Y whose rows are
%   observations and whose columns are variables and an input T-by-1 vector
%   w of weights for the observations. This function may be a valid
%   alternative to COV if observations are not all equally relevant
%   and need to be weighted according to some theoretical hypothesis or
%   knowledge.
%
%   C = WEIGHTEDCOV(Y, w) returns a positive semidefinite matrix C, i.e. all its
%   eigenvalues are non-negative.
%
%   If w = ones(size(Y, 1), 1), no difference exists between
%   WEIGHTEDCOV(Y, w) and COV(Y, 1).
%
%   REFERENCE: mathematical formulas in matrix notation are available in
%   F. Pozzi, T. Di Matteo, T. Aste,
%   "Exponential smoothing weighted correlations",
%   The European Physical Journal B, Volume 85, Issue 6, 2012.
%   DOI:10.1140/epjb/e2012-20697-x. 
%
% % ======================================================================
% % EXAMPLE
% % ======================================================================
%
% % GENERATE CORRELATED STOCHASTIC PROCESSES
%   T = 100;                                                                      % number of observations
%   N = 500;                                                                      % number of variables
%   Y = randn(T, N);                                                              % shocks from standardized normal distribution
%   Y = cumsum(Y);                                                                % correlated stochastic processes
%
% % CHOOSE EXPONENTIAL WEIGHTS
%   alpha = 2 / T;
%   w0 = 1 / sum(exp(((1:T) - T) * alpha));
%   w = w0 * exp(((1:T) - T) * alpha);                                            % weights: exponential decay
%
% % COMPUTE WEIGHTED COVARIANCE MATRIX
%   c = weightedcov(Y, w);                                                        % Weighted Covariance Matrix
%
% % ======================================================================
%
%   See also CORRCOEF, COV, STD, MEAN.
%   Check also WEIGHTEDCORRS (FE 20846) and KENDALLTAU (FE 27361)
%
% % ======================================================================
%
%   Author: Francesco Pozzi
%   E-mail: francesco.pozzi@anu.edu.au
%   Date: 15 June 2012
%
% % ======================================================================
%

% Check input
ctrl = isvector(w) & isreal(w) & ~any(isnan(w)) & ~any(isinf(w));
if ctrl
  w = w(:) / sum(w);                                                              % w is column vector
else
  error('Check w: it needs be a vector of real positive numbers with no infinite or nan values!')
end
ctrl = isreal(Y) & ~any(isnan(Y)) & ~any(isinf(Y)) & (size(size(Y), 2) == 2);
if ~ctrl
  error('Check Y: it needs be a 2D matrix of real numbers with no infinite or nan values!')
end
ctrl = length(w) == size(Y, 1);
if ~ctrl
  error('size(Y, 1) has to be equal to length(w)!')
end

[T, N] = size(Y);                                                                 % T: number of observations; N: number of variables
C = Y - repmat(w' * Y, T, 1);                                                     % Remove mean (which is, also, weighted)
C = C' * (C .* repmat(w, 1, N));                                                  % Weighted Covariance Matrix
C = 0.5 * (C + C');                                                               % Must be exactly symmetric


function f = kl(q,p)
f = q.*log(q./p) + (1-q).*log((1-q)./(1-p));


function f = normkl(mu0, S0, mu1, S1)

Si = inv(S1);

f = (trace(Si*S0) + (mu1 - mu0)'*Si*(mu1 - mu0) - log(det(S0)/det(S1)) ...
     - length(mu0))/2;


function f = entropy(mu, S)

k=length(mu);
f = k/2*(1+log(2*pi)) + log(det(S))/2;
