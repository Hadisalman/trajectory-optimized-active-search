% true model
gp_model = struct('inf',@infExact, 'mean', @meanZero, 'cov', @covSEiso, 'lik', @likGauss);
gp_para = struct('mean', [], 'cov', [0;0], 'lik', log(.1));
sn = 0.01; ell =30; sf = sqrt(1);Ncg=30;
gp_para.lik = log(sn); gp_para.cov = log([ell; sf]);


DomainBounds.xmin = 0.0;
DomainBounds.xmax = 150.0;
DomainBounds.ymin = 0.0;
DomainBounds.ymax = 150.0;
Lx = DomainBounds.xmax - DomainBounds.xmin;
Ly = DomainBounds.ymax - DomainBounds.ymin;
%Discretizing the domain
% xdel=1;
% ydel=1;
% xr=0:xdel:Lx-xdel;
% yr=0:ydel:Ly-ydel;

[X,Y]  = meshgrid(linspace(0,150,50), linspace(0,150,50));
x_gnd    = [X(:), Y(:)];
K        = gp_model.cov(gp_para.cov, x_gnd);
n        = size(x_gnd, 1);
Z = ones(size(X));
% aas params
level    = 1;
side     = 1;
highprob = .8;

% lse params
beta_t   = 9;
eps_band = .1;

num_runs = 100;
queryLen = 60;

res = 10;
[region_grid_x, region_grid_y] = meshgrid(linspace(DomainBounds.xmin,DomainBounds.xmax,res), linspace(DomainBounds.ymin,DomainBounds.ymax,res));
region_grid_x = region_grid_x(:,1:end-1);
region_grid_y = region_grid_y(1:end-1,:);
regions = [region_grid_x(:), region_grid_x(:)+Lx/res, region_grid_y(:), region_grid_y(:)+Ly/res];

recall_aas = nan(num_runs, queryLen);
recall_lse = nan(num_runs, queryLen);
recall_ei  = nan(num_runs, queryLen);
recall_unc = nan(num_runs, queryLen);
recall_rnd = nan(num_runs, queryLen);

for run_id = 1:num_runs
%    figure(1)
   sn = 0.01; ell =30; sf = sqrt(1);Ncg=30; in_noise = 0;% try 50 to notice differece 
   gp_para.lik = log(sn); gp_para.cov = [log([ell; sf])];%;in_noise];

    sGT = chol(K + exp(2*gp_para.lik) * eye(n))' * randn(n,1);
    sGT=max(sGT,0);
    sGT=5*sGT/(max(sGT));
    y_gnd = sGT;
    surface(X,Y,Z,reshape(sGT,size(X)))
    shading interp
    colorbar
%     
    sn = 0.01; ell =15; sf = sqrt(1);Ncg=30; in_noise = 0;% try 50 to notice differece 
    gp_para.lik = log(sn); gp_para.cov = [log([ell; sf])];%;in_noise];

%    y_gnd = chol(K + exp(2*gp_para.lik) * eye(n))' * randn(n,1);
% 
%    surface(X1,X2,0*X2,reshape(y_gnd,50,50))
%    shading interp
%    colorbar
%    drawnow
%    ground truth
   gnd = ActiveAreaSearch(gp_model, gp_para, x_gnd, regions, level, side, highprob);
   region_outcome_gnd = gnd.update(x_gnd, y_gnd);


   % initialize search algorithms
   aas = ActiveAreaSearch        (gp_model, gp_para, x_gnd, regions, level, side, highprob);
   lse = ActiveLevelSetEstimation(gp_model, gp_para, x_gnd, level, beta_t, eps_band);
   unc = UncertaintySampling     (gp_model, gp_para, x_gnd);
   ei  = ExpectedImprovement     (gp_model, gp_para, x_gnd);
   rnd = RandomSampling          (gp_model, gp_para, x_gnd);

   % initialize measures
   region_measure_lse = ActiveAreaSearch(gp_model, gp_para, x_gnd, regions,level, side, highprob);
   region_measure_unc = ActiveAreaSearch(gp_model, gp_para, x_gnd, regions, level, side, highprob);
   region_measure_ei  = ActiveAreaSearch(gp_model, gp_para, x_gnd, regions, level, side, highprob);
   region_measure_rnd = ActiveAreaSearch(gp_model, gp_para, x_gnd, regions, level, side, highprob);

   for query_count = 0:queryLen-1
      % aas
      u = aas.utility();
      [~, ind] = max_tiebreak(u,[],false);
      aas.update(x_gnd(ind, :), y_gnd(ind, :));
      
      [~, ~, f, fs2] = aas.predict_points(x_gnd);
%         figure(2)
%       surface(X1,X2,0*X2,reshape(f,50,50))
%       shading interp
%       colorbar
%       drawnow
%       pause()
      
      recall_aas(run_id, query_count+1) = (0+aas.cumfound>0)'*region_outcome_gnd / sum(region_outcome_gnd);

      % lse
      u = lse.utility();
      [~, ind] = max_tiebreak(u,[],false);
      lse.update(x_gnd(ind, :), y_gnd(ind, :));

      region_measure_lse.update(x_gnd(ind, :), y_gnd(ind, :));

      recall_lse(run_id, query_count+1) = (0+region_measure_lse.cumfound>0)'*region_outcome_gnd  / sum(region_outcome_gnd);

      % unc
      u = unc.utility();
      [~, ind] = max_tiebreak(u,[],false);
      unc.update(x_gnd(ind, :), y_gnd(ind, :));

      region_measure_unc.update(x_gnd(ind, :), y_gnd(ind, :));

      recall_unc(run_id, query_count+1) = (0+region_measure_unc.cumfound>0)'*region_outcome_gnd / sum(region_outcome_gnd);

      % ei
      u = ei.utility();
      [~, ind] = max_tiebreak(u,[],false);
      ei.update(x_gnd(ind, :), y_gnd(ind, :));
      [~, ~, f, fs2] = ei.predict_points(x_gnd);
 
%       figure(2)
%       surface(X1,X2,0*X2,reshape(f,50,50))
%       shading interp
%       colorbar
%       drawnow

      region_measure_ei.update(x_gnd(ind, :), y_gnd(ind, :));

      recall_ei(run_id, query_count+1) = (0+region_measure_ei.cumfound>0)'*region_outcome_gnd / sum(region_outcome_gnd);

      % rand
      u = rnd.utility();
      [~, ind] = max_tiebreak(u,[],false);
      rnd.update(x_gnd(ind, :), y_gnd(ind, :));

      region_measure_rnd.update(x_gnd(ind, :), y_gnd(ind, :));

      recall_rnd(run_id, query_count+1) = (0+region_measure_rnd.cumfound>0)'*region_outcome_gnd / sum(region_outcome_gnd);
   end

   [recall_aas(run_id, end), recall_lse(run_id, end), recall_unc(run_id, end), recall_ei(run_id, end), recall_rnd(run_id, end)]

end
%%
figure;

% errorbar(1:queryLen, mean(recall_aas, 1), zeros(size(mean(recall_aas, 1))))
% hold on
% errorbar(1:queryLen, mean(recall_lse, 1), zeros(size(mean(recall_aas, 1))))
% errorbar(1:queryLen, mean(recall_unc, 1), zeros(size(mean(recall_aas, 1))))
% errorbar(1:queryLen, mean(recall_ei, 1), zeros(size(mean(recall_aas, 1))))
% errorbar(1:queryLen, mean(recall_rnd, 1), zeros(size(mean(recall_aas, 1))))
% legend('aas','lse','unc','ei','rand');
% title('recall curves (mean and standard error from 100 runs)');
% grid on;


errorbar(1:queryLen, mean(recall_aas, 1), std(recall_aas, 1)/sqrt(num_runs))
hold on
errorbar(1:queryLen, mean(recall_lse, 1), std(recall_lse, 1)/sqrt(num_runs))
errorbar(1:queryLen, mean(recall_unc, 1), std(recall_unc, 1)/sqrt(num_runs))
errorbar(1:queryLen, mean(recall_ei, 1), std(recall_ei, 1)/sqrt(num_runs))
errorbar(1:queryLen, mean(recall_rnd, 1), std(recall_rnd, 1)/sqrt(num_runs))
legend('aas','lse','unc','ei','rand');
title('recall curves (mean and standard error from 100 runs)');
grid on;

