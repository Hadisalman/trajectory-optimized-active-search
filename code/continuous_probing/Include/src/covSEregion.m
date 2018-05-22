function [omega, Z] = covSEregion(hyp, regions, x)
% require regions = [xmin, xmax, ymin, ymax]

loc_dim = size(regions, 2)/2;

sf2     = exp(2 * hyp(end)); 
ell     = exp(    hyp(1:end-1));
if isscalar(ell)
  % convert covSEiso to covSEard
  ell = repmat(ell, 1, loc_dim);
end

% necessary normalizers
logZ   = sum ( ...
  .5*log(2*pi) * ones(1,loc_dim) + log( ell) ...
  );
gamma2 = exp( log( sf2 ) + logZ );

[n, ~] = size(x);

% initialize
omega = nan(size(regions, 1), n);
Z     = nan(size(regions, 1), 1);


for groupid = 1:size(regions, 1)
  
  %----------------- parameters ------------------  
  lims_lower = regions(groupid, 1:2:size(regions, 2));
  lims_upper = regions(groupid, 2:2:size(regions, 2));

  assert(all(lims_lower<lims_upper));
  
  area = prod(lims_upper - lims_lower);
  
  %----------------- omega ------------------
  if n>0
    omega(groupid, :) = gamma2/area * prod( ...
      normcdf( bsxfun(@rdivide, bsxfun(@minus, lims_upper(:)', x), ell(:)' ) ) ...
      - normcdf( bsxfun(@rdivide, bsxfun(@minus, lims_lower(:)', x), ell(:)') ) , ...
      2)';
  end

  %----------------- Z if needed ------------------
  if nargout >= 2 
    
    integrals = nan(size(ell));
    
    for d=1:loc_dim
      integrand = @(x) normcdf((lims_upper(d) - x)./ell(d)) ...
        - normcdf((lims_lower(d) - x)./ell(d)) ;
      integrals(d) = quad(integrand, lims_lower(d), lims_upper(d));
    end
    
    Z(groupid) = gamma2/area^2 * prod( integrals );
        
  end

end
