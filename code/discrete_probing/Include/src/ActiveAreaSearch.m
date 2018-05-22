classdef ActiveAreaSearch < ActiveGP
  properties
  	regions     % g * 2d
  	level       % scalar
  	side        % scalar
  	highprob    % scalar
    cumfound    % g * 1
    utility_omit_positive_regions = true % only compute expected reward on undiscovered regions
    pool_locs   % ns * d
  end

  properties (Transient=true)
  	Z           % g * 1
  	omega       % g * n
    alpha       % g * 1, region posterior mean
    beta2       % g * 1, region posterior variance
    region_gpml_alpha   % g * 1, for expected reward computation
    omegas
    tV
  end

  methods
    function self = ActiveAreaSearch(gp_model, gp_para, pool_locs, regions, level, side, highprob)

      self = self@ActiveGP(gp_model, gp_para);

      self.regions   = regions;
      self.level     = level;
      self.side      = side;
      self.highprob  = highprob;
      self.cumfound  = zeros(size(self.regions, 1), 1);

      [~, self.Z]    = covSEregion(self.gp_para.cov, self.regions, []);
      self.alpha     = 0 * self.Z;
      self.beta2     = self.Z;

      self.pool_locs = pool_locs;
    end

    function set.pool_locs(self, pool_locs)
      [omegas, tV] = self.pool_locs_compute_stats(pool_locs);
      self.pool_locs = pool_locs;
      self.omegas = omegas;
      self.tV = tV;
    end

    function [new_found, Tg] = update(self, new_locs, new_vals)
      update@ActiveGP(self, new_locs, new_vals);
      [new_found, Tg] = self.update_regions();
    end

    function [new_found, Tg] = update_regions(self)
      assert(~isempty(self.collected_locs))

      % update region kernel cross terms
      self.omega = [self.omega, covSEregion(self.gp_para.cov, self.regions, ...
        self.collected_locs((size(self.omega, 2)+1):end, :))];

      % compute and cache region stats
      self.alpha = self.omega * self.gp_post.alpha;

      V  = self.R' \ self.omega';
      self.beta2 = self.Z - sum(V.*V,1)';     % predictive variances

      self.region_gpml_alpha = solve_chol(self.R, self.omega');

      assert(all(self.beta2>0));

      % compute region outcomes using cached region stats
      [new_found, Tg] = self.region_rewards();
      self.cumfound = self.cumfound + new_found;
    end

    function [new_found, Tg] = region_rewards(self)
      Tg = normcdf(self.side .* (self.alpha - self.level) ./ sqrt(self.beta2));
      new_found = 0 + (Tg > self.highprob);
    end 

    function [omegas, tV] = pool_locs_compute_stats(self, pool_locs)
      omegas = covSEregion(self.gp_para.cov, self.regions, pool_locs);
      [~, tV] = self.predict_points(pool_locs);
    end

    function [u, ug] = utility(self, pool_locs)
      if nargin<2
        pool_locs = self.pool_locs;
        omegas    = self.omegas;
        tV        = self.tV;
      else
        [omegas, tV] = self.pool_locs_compute_stats(pool_locs);
      end

      ns = size(pool_locs, 1);
      g  = size(self.regions, 1);

      if ~isempty(self.collected_locs)
        ks = self.gp_model.cov(self.gp_para.cov, self.collected_locs, pool_locs);
        omega_Vinv_ks = self.region_gpml_alpha' * ks;
      else
        % avoid calling cov on an empty matrices without the correct shape
        omega_Vinv_ks = 0;
      end
      
      tnu = (abs(omegas - omega_Vinv_ks) ./ repmat(sqrt(tV)', g, 1));  % g * n
      tbeta2 = repmat(self.beta2, 1, ns) - tnu.^2;  % g * n

      assert(all(all(tbeta2>0)));

      alpham = repmat(self.alpha, 1, ns);

      ug = normcdf( (self.side .* (alpham - self.level) - sqrt(tbeta2).*norminv(self.highprob)) ./ tnu );

      if self.utility_omit_positive_regions
        u = sum(ug .* repmat(self.cumfound==0, 1, ns), 1);
      else
        u = sum(ug, 1);
      end
    end
  end
end
