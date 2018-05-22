classdef ExpectedImprovement < ActiveGP
    properties
        best_observed
        best_location
        pool_locs
    end
    
    methods
        function self = ExpectedImprovement(gp_model, gp_para, pool_locs)
            self = self@ActiveGP(gp_model, gp_para);
            self.best_observed = -inf;
            if nargin>=3
                self.pool_locs = pool_locs;
            end
        end
        
        function [best_observed] = update(self, locations, values)
            update@ActiveGP(self, locations, values);
            
            if max(values) > self.best_observed
                [~, ind] = max(values);
                self.best_observed = max(values);
                self.best_location = locations(ind, :);
            end
        end
        
        function [EI] = utility(self, pool_locs)
            if nargin<2
                pool_locs = self.pool_locs;
            end
            [~, ~, f, fs2] = self.predict_points(pool_locs);
            eps=0.01;            
            sigma = sqrt(fs2);
            Z = (f - self.best_observed-eps) ./ sigma;
            EI = (f - self.best_observed) .* normcdf(Z) + sigma .* normpdf(Z);
        end
    end
end
