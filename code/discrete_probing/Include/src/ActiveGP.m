classdef ActiveGP < handle
  properties
   gp_model
   gp_para
   collected_locs  % n * d
   collected_vals  % n * 1
  end

  properties (Transient=true)
    gp_post
    R
  end

  methods
    function self = ActiveGP(gp_model, gp_para)

      assert( exp(2*gp_para.lik) > 1e-6 )  % gpml uses different representations otherwise
      assert( strcmp(func2str(gp_model.mean), 'meanZero') ) % simplify model

      self.gp_model   = gp_model;
      
      if ~isfield(gp_para, 'mean')
        gp_para.mean  = [];
      end
      self.gp_para    = gp_para;
    end

    function [ymu, ys2, fmu, fs2] = predict_points(self, new_x)
      if isempty(self.collected_locs)
        fmu = self.gp_model.mean(self.gp_para.mean, new_x);
        fs2 = self.gp_model.cov (self.gp_para.cov , new_x, 'diag');
        ymu = fmu;
        ys2 = fs2 + exp(2*self.gp_para.lik);
      else
        [ymu, ys2, fmu, fs2] = gp(self.gp_para, self.gp_model.inf, self.gp_model.mean, self.gp_model.cov, self.gp_model.lik, self.collected_locs, self.gp_post, new_x);
      end
    end

    function [] = update(self, locations, values)
      if isempty(self.collected_locs)
        [~, ~, self.gp_post] = gp(self.gp_para, self.gp_model.inf, self.gp_model.mean, self.gp_model.cov, self.gp_model.lik, locations, values(:));
        self.R = self.gp_post.L .* exp(self.gp_para.lik);
      else
        self.gp_post = update_posterior(self.gp_para, self.gp_model.mean, {self.gp_model.cov}, self.collected_locs, self.gp_post, locations, values(:));
        self.R = self.gp_post.L .* exp(self.gp_para.lik);
      end

      self.collected_locs = [self.collected_locs; locations ];
      self.collected_vals = [self.collected_vals; values(:) ];
    end
  end
end
