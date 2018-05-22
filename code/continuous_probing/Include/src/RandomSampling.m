classdef RandomSampling < ActiveGP
  properties
    pool_locs
  end
  
  methods
    function self = RandomSampling(gp_model, gp_para, pool_locs)
      self = self@ActiveGP(gp_model, gp_para);
      if nargin>=3
        self.pool_locs = pool_locs;
      end
    end

    function [u] = utility(self, pool_locs)
      if nargin<2
         pool_locs = self.pool_locs;
      end
      u = rand(size(pool_locs, 1), 1);
    end
  end
end
