%Authors: Hadi Salman, Elif Ayvali, and Yifei Ma
%Please acknowledge:
%GPML toolkit by Carl Edward Rasmussen and Hannes Nickisch
%GPML_EXTENSIONS by Roman Garnett
%Cite:
%Active Area Search paper by Yifei Ma et al
%Active Level Sets Alkis Gotovos et al

addpath(genpath(pwd))
clc;close all;clear
load 2tumors_experiment.mat
% load 1tumor_1artery.mat

x1 = reshape(x_gnd(:,1), [51,51])';
x2 = reshape(x_gnd(:,2), [51,51])';
x3 = reshape(x_gnd(:,3), [51,51])';
y = reshape(y_gnd, [51,51])';
x_gnd = [x1(:), x2(:), x3(:)];
y_gnd = y(:);

surface(x1,x2,x3,y)
title('ground truth stiffness map')
shading interp
axis equal tight
colorbar

L1 = x1(1,end) - x1(1,1);
L2 = x2(end,1) - x2(1,1);
[region_grid_x, region_grid_y] = meshgrid(linspace(0,L1,11), linspace(0,L2,11));
region_grid_x = region_grid_x(:,1:end-1);
region_grid_y = region_grid_y(1:end-1,:);
regions = [region_grid_x(:), region_grid_x(:)+L1/11, region_grid_y(:), region_grid_y(:)+L2/11];

x_gnd=x_gnd(:,1:2);
gp_model = struct('inf',@infExact, 'mean', @meanZero, 'cov', @covSEiso, 'lik', @likGauss);
sn = 0.01; ell =8; sf = sqrt(1);Ncg=30;
gp_para.lik = log(sn); gp_para.cov = log([ell; sf]);
%%

queryLen = 60;
num_runs = 10;
recall_aas = nan(num_runs, queryLen);
recall_lse = nan(num_runs, queryLen);
recall_ei  = nan(num_runs, queryLen);
recall_unc = nan(num_runs, queryLen);
recall_rnd = nan(num_runs, queryLen);

x_shape  = [51,51];%:size(x_gnd)
side     = 1;
beta_t   = 9;
eps_band = .1;
level    = 1;

gnd = ActiveAreaSearch(gp_model, gp_para, x_gnd, regions, level, side, highprob);
found = gnd.update(x_gnd, y_gnd);
region_outcome_gnd = found;
%     figure(1); clf;
%     [~, ~, f, fs2] = gnd.predict_points(x_gnd);
%     plot_demo(x_shape, x_gnd, f, fs2, [], regions, level, gnd.alpha, gnd.beta2, found, 'plotTailGaussian', false);
%     color_line3(x_gnd(:,1),x_gnd(:,2),ones(length(x_gnd(:,2)),1),y_gnd,'.')
%     set(gcf, 'Position', [000, 100, 280, 280]);axis equal;
%     axis equal
%     title('ground truth');
%     drawnow;
for run_id = 1: num_runs
    display(run_id)

    %%
    % ----------------------- aas ----------------------------
    aas = ActiveAreaSearch(gp_model, gp_para, x_gnd, regions, level, side, highprob);
%     figure(2); clf;hold on;
%     set(gcf, 'Position', [300, 100, 280, 280]); axis equal;
%     title('active area search');
    for query_count = 0:queryLen-1
        u = aas.utility();
        if(query_count ==0)
            ind = randi(size(x_gnd,1),1,1);
            ind_start=ind;
        else
        [~, ind] = max_tiebreak(u);
        end
        found = aas.update(x_gnd(ind, :), y_gnd(ind, :));
        recall_aas(run_id, query_count+1) = (0+aas.cumfound>0)'*region_outcome_gnd / sum(region_outcome_gnd);
        %to see results incrementally:
%         [~, ~, f, fs2] = aas.predict_points(x_gnd);
%         plot_demo(x_shape, x_gnd, f, fs2, aas.collected_locs, regions, level, aas.alpha, aas.beta2, found, 'plotTailGaussian', false);
%         set(gcf, 'Position', [300, 100, 280, 280]); axis equal;
%         title('active area search');
    %     drawnow;
%         pause(0.1);

    end
    %plot all
%     figure(2); clf;hold on;
%     [~, ~, f, fs2] = aas.predict_points(x_gnd);
%     plot_demo(x_shape, x_gnd, f, fs2, aas.collected_locs, regions, level, aas.alpha, aas.beta2, found, 'plotTailGaussian', false);
    % color_line3(x_gnd(:,1),x_gnd(:,2),ones(length(x_gnd(:,1)),1),f,'.');hold on
    % scatter(aas.collected_locs(:,1), aas.collected_locs(:,2), 20, 'mo', 'filled');
%     set(gcf, 'Position', [300, 100, 280, 280]); axis equal;
%     title('active area search');
%     drawnow;
    %%
    % ----------------------- lse ----------------------------
    lse = ActiveLevelSetEstimation(gp_model, gp_para, x_gnd, level, beta_t, eps_band);
%     figure(3); clf;hold on;
%     set(gcf, 'Position', [600, 100, 280, 280]);axis equal;
%     title('active level set estimation');
    for query_count = 0:queryLen-1
        u = lse.utility();
        if query_count==0 %for all algorithms to start from same point
            ind = ind_start;
        else
        [~, ind] = max_tiebreak(u);
        end
        lse.update(x_gnd(ind, :), y_gnd(ind, :));

    %    %to see results incrementally:
        region_measure_lse = ActiveAreaSearch(gp_model, gp_para, x_gnd, regions, level, side, highprob);
        found = region_measure_lse.update(lse.collected_locs, lse.collected_vals);
        recall_lse(run_id, query_count+1) = (0+region_measure_lse.cumfound>0)'*region_outcome_gnd  / sum(region_outcome_gnd);
%         [~, ~, f, fs2] = region_measure_lse.predict_points(x_gnd);
%         plot_demo(x_shape, x_gnd, f, fs2, lse.collected_locs, regions, level, region_measure_lse.alpha, region_measure_lse.beta2, found, 'plotTailGaussian', false);
    %     drawnow;
%         pause(0.1);
    end

    %plot all results
%     region_measure_lse = ActiveAreaSearch(gp_model, gp_para, x_gnd, regions, level, side, highprob);
%     found = region_measure_lse.update(lse.collected_locs, lse.collected_vals);

%     figure(3); clf;hold on;
%     [~, ~, f, fs2] = region_measure_lse.predict_points(x_gnd);
%     plot_demo(x_shape, x_gnd, f, fs2, lse.collected_locs, regions, level, region_measure_lse.alpha, region_measure_lse.beta2, found, 'plotTailGaussian', false);
%     % color_line3(x_gnd(:,1),x_gnd(:,2),ones(length(x_gnd(:,1)),1),f,'.');hold on
%     % scatter(lse.collected_locs(:,1), lse.collected_locs(:,2), 20, 'mo', 'filled');
%     set(gcf, 'Position', [600, 100, 280, 280]);axis equal;
%     title('active level set estimation');
%     drawnow;

    %%
    % ----------------------- unc ----------------------------
    unc = UncertaintySampling(gp_model, gp_para, x_gnd);
%     figure(4); clf;hold on;
%     set(gcf, 'Position', [000, 400, 280, 280]);axis equal;
%     title('uncertainty sampling');

    for query_count = 0:queryLen-1
        u = unc.utility();
        if query_count==0 %for all algorithms to start from same point
            ind = ind_start;
        else
        [~, ind] = max_tiebreak(u);
        end
        unc.update(x_gnd(ind, :), y_gnd(ind, :));

        %to see results incrementally:
        region_measure_unc = ActiveAreaSearch(gp_model, gp_para, x_gnd, regions, level, side, highprob);
        found = region_measure_unc.update(unc.collected_locs, unc.collected_vals);
        recall_unc(run_id, query_count+1) = (0+region_measure_unc.cumfound>0)'*region_outcome_gnd / sum(region_outcome_gnd);
%         [~, ~, f, fs2] = region_measure_unc.predict_points(x_gnd);
%         plot_demo(x_shape, x_gnd, f, fs2, unc.collected_locs, regions, level, region_measure_unc.alpha, region_measure_unc.beta2, found, 'plotTailGaussian', false);
    %     drawnow
%         pause(0.1);

    end

%     region_measure_unc = ActiveAreaSearch(gp_model, gp_para, x_gnd, regions, level, side, highprob);
%     found = region_measure_unc.update(unc.collected_locs, unc.collected_vals);

%     figure(4); clf;hold on;
%     [~, ~, f, fs2] = region_measure_unc.predict_points(x_gnd);
%     plot_demo(x_shape, x_gnd, f, fs2, unc.collected_locs, regions, level, region_measure_unc.alpha, region_measure_unc.beta2, found, 'plotTailGaussian', false);
%     % color_line3(x_gnd(:,1),x_gnd(:,2),ones(length(x_gnd(:,1)),1),f,'.')
%     % scatter(unc.collected_locs(:,1), unc.collected_locs(:,2), 20, 'mo', 'filled');
% 
%     set(gcf, 'Position', [000, 400, 280, 280]);axis equal;
%     title('uncertainty sampling');
%     drawnow;
    %%

    % ----------------------- ei ----------------------------
    ei = ExpectedImprovement(gp_model, gp_para, x_gnd);
%     figure(5); clf;hold on;
%     set(gcf, 'Position', [300, 400, 280 , 280]);axis equal;
%     title('expected improvement');
    for query_count = 0:queryLen-1
        u = ei.utility();
        if query_count==0 %for all algorithms to start from same point
            ind = ind_start;
        else
        [~, ind] = max_tiebreak(u);
        end
        ei.update(x_gnd(ind, :), y_gnd(ind, :));
        %to see results incrementally:
        region_measure_ei = ActiveAreaSearch(gp_model, gp_para, x_gnd, regions, level, side, highprob);
        found = region_measure_ei.update(ei.collected_locs, ei.collected_vals);
        recall_ei(run_id, query_count+1) = (0+region_measure_ei.cumfound>0)'*region_outcome_gnd / sum(region_outcome_gnd);

%         [~, ~, f, fs2] = region_measure_ei.predict_points(x_gnd);
%         plot_demo(x_shape, x_gnd, f, fs2, ei.collected_locs, regions, level, region_measure_ei.alpha, region_measure_ei.beta2, found, 'plotTailGaussian', false);
    %     drawnow;
%         pause(0.1)

    end

%     region_measure_ei = ActiveAreaSearch(gp_model, gp_para, x_gnd, regions, level, side, highprob);
%     found = region_measure_ei.update(ei.collected_locs, ei.collected_vals);

%     figure(5); clf;hold on;
%     [~, ~, f, fs2] = region_measure_ei.predict_points(x_gnd);
%     plot_demo(x_shape, x_gnd, f, fs2, ei.collected_locs, regions, level, region_measure_ei.alpha, region_measure_ei.beta2, found, 'plotTailGaussian', false);
    % color_line3(x_gnd(:,1),x_gnd(:,2),ones(length(x_gnd(:,1)),1),f,'.')
    % scatter(ei.collected_locs(:,1), ei.collected_locs(:,2), 20, 'mo', 'filled');
%     set(gcf, 'Position', [300, 400, 280 , 280]);axis equal;
%     title('expected improvement');
%     drawnow;

    %%
    % ----------------------- rnd ----------------------------
    rnd = RandomSampling(gp_model, gp_para, x_gnd);
%     figure(6); clf;hold on;
%     set(gcf, 'Position', [600, 400, 280, 280]);axis equal;
%     title('random sampling');

    for query_count = 0:queryLen-1
        u = rnd.utility();
        if query_count==0 %for all algorithms to start from same point
            ind = ind_start;
        else
        [~, ind] = max_tiebreak(u);
        end
        rnd.update(x_gnd(ind, :), y_gnd(ind, :));
        %to see results incrementally:
        region_measure_rnd = ActiveAreaSearch(gp_model, gp_para, x_gnd, regions, level, side, highprob);
        found = region_measure_rnd.update(rnd.collected_locs, rnd.collected_vals);
        recall_rnd(run_id, query_count+1) = (0+region_measure_rnd.cumfound>0)'*region_outcome_gnd / sum(region_outcome_gnd);

%         [~, ~, f, fs2] = region_measure_rnd.predict_points(x_gnd);
%         plot_demo(x_shape, x_gnd, f, fs2, rnd.collected_locs, regions, level, region_measure_rnd.alpha, region_measure_rnd.beta2, found, 'plotTailGaussian', false);
%     %     drawnow;
%         pause(0.1);
    end

%     region_measure_rnd = ActiveAreaSearch(gp_model, gp_para, x_gnd, regions, level, side, highprob);
%     found = region_measure_rnd.update(rnd.collected_locs, rnd.collected_vals);
% 
%     figure(6); clf;hold on;
%     [~, ~, f, fs2] = region_measure_rnd.predict_points(x_gnd);
%     plot_demo(x_shape, x_gnd, f, fs2, rnd.collected_locs, regions, level, region_measure_rnd.alpha, region_measure_rnd.beta2, found, 'plotTailGaussian', false);
%     % color_line3(x_gnd(:,1),x_gnd(:,2),ones(length(x_gnd(:,1)),1),f,'.')
%     % scatter(rnd.collected_locs(:,1), rnd.collected_locs(:,2), 20, 'mo', 'filled');
%     set(gcf, 'Position', [600, 400, 280, 280]);axis equal;
%     title('random sampling');
%     drawnow;

end

%% Plotting recall measure

% zeros(size(mean(recall_aas, 1)))
% figure;
% 
% errorbar(1:queryLen, mean(recall_aas, 1), zeros(size(mean(recall_aas, 1))))
% hold on
% errorbar(1:queryLen, mean(recall_lse, 1), zeros(size(mean(recall_aas, 1))))
% errorbar(1:queryLen, mean(recall_unc, 1), zeros(size(mean(recall_aas, 1))))
% errorbar(1:queryLen, mean(recall_ei, 1), zeros(size(mean(recall_aas, 1))))
% errorbar(1:queryLen, mean(recall_rnd, 1), zeros(size(mean(recall_aas, 1))))
% legend('aas','lse','unc','ei','rand');
% title('recall curves (mean and standard error from 100 runs)');
% grid on;

figure;
hold on
errorbar(1:queryLen, mean(recall_aas, 1), std(recall_aas, 1)/sqrt(num_runs))
errorbar(1:queryLen, mean(recall_lse, 1), std(recall_lse, 1)/sqrt(num_runs))
errorbar(1:queryLen, mean(recall_unc, 1), std(recall_unc, 1)/sqrt(num_runs))
errorbar(1:queryLen, mean(recall_ei, 1), std(recall_ei, 1)/sqrt(num_runs))
errorbar(1:queryLen, mean(recall_rnd, 1), std(recall_rnd, 1)/sqrt(num_runs))
legend('aas','lse','unc','ei','rand');
title('recall curves (mean and standard error from 100 runs)');
grid on;
