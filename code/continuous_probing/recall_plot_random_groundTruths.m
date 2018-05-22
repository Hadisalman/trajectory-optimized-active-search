%Authors: Hadi Salman, Elif Ayvali, and Yifei Ma

%% Initialization
clc;close all;clear;
opt = [];%initialize options

opt.gp_model = struct('inf',@infExact, 'mean', @meanZero, 'cov', @covSEiso, 'lik', @likGauss);
sn = 0.01; ell =30; sf = sqrt(1);Ncg=30; in_noise = 0;% try 50 to notice differece 
opt.gp_para.lik = log(sn); opt.gp_para.cov = [log([ell; sf])];%;in_noise];
draw = 0;

% aas params
level    = 1;
side     = 1;
highprob = .8;

% lse params
beta_t   = 9;
eps_band = .1;

%%%%%%%%%%%%%%%%Definition of Variables%%%%%%%%%%%%%%%%%%%%
%xss (m x d):3-d organ grid
%sGT (m x 1):stiffness ground truth
%AF (m x 1) :acquisition function
%posGP (1 x d):probed points
%kGP (n x 1):stiffness at the probed points
%posGP_next (1 x d):points to be probed along the trajectory
%ymu (m x 1) :GP mean
%ys2 (m x 1) :GP uncertainty
%xs_traj ((d+1) x L): trajectory with L ( d-dimensional position and orientation)
%Yifei: only main.m and EvaluateTrajUtility.m has GP related stuff
%--------------------------------------------------------%
%%%%%%%%%%%%%%%%Setting domain bounds%%%%%%%%%%%%%%%%%%%%%
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
K        = opt.gp_model.cov(opt.gp_para.cov, x_gnd);
n        = size(x_gnd, 1);

num_runs = 100;
queryLen = 60;

recall_aas = nan(num_runs, queryLen);
recall_lse = nan(num_runs, queryLen);
recall_ei  = nan(num_runs, queryLen);
recall_unc = nan(num_runs, queryLen);
recall_rnd = nan(num_runs, queryLen);

% [region_X1, region_X2] = meshgrid(0:9, 0:9);
% regions  = [region_X1(:), region_X1(:)+1, region_X2(:), region_X2(:)+1];
res = 10;
[region_grid_x, region_grid_y] = meshgrid(linspace(DomainBounds.xmin,DomainBounds.xmax,res), linspace(DomainBounds.ymin,DomainBounds.ymax,res));
region_grid_x = region_grid_x(:,1:end-1);
region_grid_y = region_grid_y(1:end-1,:);
regions = [region_grid_x(:), region_grid_x(:)+Lx/res, region_grid_y(:), region_grid_y(:)+Ly/res];

%--------------------------------------------------------%
%%%%%%%%%%%%%%%%Generate Stiffness Map%%%%%%%%%%%%%%%%%%%%%
% addnoise=0;
% [X,Y,sGT] = GenerateStiffnessMap(xr,yr,addnoise);
% sGT = 4*reshape(sGT,size(X));
Z   = ones(size(X)); %for 3D problem
xss=[X(:),Y(:),Z(:)];
%(3d point, scalar stiffness at the point)
% datafull = [xss,sGT(:)];
%%
for run_id = 1:num_runs

   sn = 0.01; ell =30; sf = sqrt(1);Ncg=30; in_noise = 0;% try 50 to notice differece 
   opt.gp_para.lik = log(sn); opt.gp_para.cov = [log([ell; sf])];%;in_noise];
    if run_id==1
    sGT = chol(K + exp(2*opt.gp_para.lik) * eye(n))' * randn(n,1);
    sGT=max(sGT,0);
    sGT=5*sGT/(max(sGT));
    y_gnd = sGT;
    end
%     surface(X,Y,Z,reshape(sGT,size(X)))
%     shading interp
%     colorbar
%     
    sn = 0.01; ell =15; sf = sqrt(1);Ncg=30; in_noise = 0;% try 50 to notice differece 
    opt.gp_para.lik = log(sn); opt.gp_para.cov = [log([ell; sf])];%;in_noise];

    gnd = ActiveAreaSearch(opt.gp_model, opt.gp_para, x_gnd, regions, level, side, highprob);
   region_outcome_gnd = gnd.update(x_gnd, y_gnd);


   % initialize search algorithms
   aas = ActiveAreaSearch        (opt.gp_model, opt.gp_para, x_gnd, regions, level, side, highprob);
   lse = ActiveLevelSetEstimation(opt.gp_model, opt.gp_para, x_gnd, level, beta_t, eps_band);
   unc = UncertaintySampling     (opt.gp_model, opt.gp_para, x_gnd);
   ei  = ExpectedImprovement     (opt.gp_model, opt.gp_para, x_gnd);
   rnd = RandomSampling          (opt.gp_model, opt.gp_para, x_gnd);

   % initialize measures
   region_measure_lse = ActiveAreaSearch(opt.gp_model, opt.gp_para, x_gnd, regions,level, side, highprob);
   region_measure_unc = ActiveAreaSearch(opt.gp_model, opt.gp_para, x_gnd, regions, level, side, highprob);
   region_measure_ei  = ActiveAreaSearch(opt.gp_model, opt.gp_para, x_gnd, regions, level, side, highprob);
   region_measure_rnd = ActiveAreaSearch(opt.gp_model, opt.gp_para, x_gnd, regions, level, side, highprob);
idxinit = randi(size(xss,1),1,1); % any random initial point
for alg = 1:4
%--------------------------------------------------------%
%%%%%%%%%%Initialize Prior Acquisition Function%%%%%%%%%%%%%
opt.gp.AF = zeros(size(X));
%--------------------------------------------------------%
%%%%%%%%%%%%%%initialize robot position%%%%%%%%%%%%%%%%%%%%
%Get te index of the position where AF is max
% [~,idxinit]=max(opt.gp.AF(:));

%--------------------------------------------------------%
%%%%%%%%%%%%%%%%initialize CE planner%%%%%%%%%%%%%%%%%%%%%%%
[opt] = initialize_gen_traj_CE(opt,DomainBounds);
opt.planner=1;
if(opt.planner==1)
    display('Planner:Utility Maximization')
end

if(opt.planner==2)
    display('Planner:Ergodic Coverage')
end

%If you add a new method, don't forget to update EvaluateUtility.m

opt.method=alg;
if(opt.method==1)
    display('Method: Expected improvement')
end

if(opt.method==2)
    display('Method:Variance Reduction')
end

if(opt.method==1)
    display('Method: Level Set Estimation')
end

if(opt.method==4)
    display('Method: Active Area Search')
end

opt.xi = [xss(idxinit,:)'; -0.0*pi];%[pos;orientation]
traj_save=opt.xi; %save trajectories
%--------------------------------------------------------%
%%%%%%%%%%%%%%%%initialize ergodicity%%%%%%%%%%%%%%%%%%%%%%%
if(opt.planner==2)
    opt.erg.mu=opt.gp.AF./sum(sum(opt.gp.AF));
    % Number of wave-numbers to be used
    Nk = 10;%%
    opt.erg.Nkx = Nk;
    opt.erg.Nky = Nk;
    %Undo Fourier reflection
%     opt.erg.mu=flipud(opt.erg.mu);
%     opt.erg.mu=imrotate(opt.erg.mu,-90);
    [opt.erg.muk] = GetFourierCoeff(opt,X,Y);
    opt.erg.Ck = zeros(Nk,Nk);
end
%--------------------------------------------------------%
%%%%%%%%%%%%%%%%initialize GP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opt.gp.posGP=[];opt.gp.kGP=[];sample_rate=5;
%GP parameters
opt.gp.posGP=opt.xi(1:opt.dim,:)';
opt.gp.kGP=EvaluateStiffnessKnn(opt.gp.posGP,xss,sGT);

if(opt.method==3)
    opt.level = 1;
    opt.beta_t = 9;
    opt.eps_band = .1;
    opt.method_obj = ActiveLevelSetEstimation(opt.gp_model, opt.gp_para, xss, opt.level, opt.beta_t, opt.eps_band);
end

if(opt.method==4)
    opt.level=1;
    res = 20;
    [region_grid_x, region_grid_y] = meshgrid(linspace(DomainBounds.xmin,DomainBounds.xmax,res), linspace(DomainBounds.ymin,DomainBounds.ymax,res));
    region_grid_x = region_grid_x(:,1:end-1);
    region_grid_y = region_grid_y(1:end-1,:);
    opt.regions = [region_grid_x(:), region_grid_x(:)+Lx/res, region_grid_y(:), region_grid_y(:)+Ly/res];

%     [region_grid_x, region_grid_y] = meshgrid(0:10:130, 0:10:130);
%     opt.regions = [region_grid_x(:), region_grid_x(:)+20, region_grid_y(:), region_grid_y(:)+20];
%     opt.regions(:,5) = -0.1;
%     opt.regions(:,6) = 0.1;
    opt.method_obj = ActiveAreaSearch(opt.gp_model, opt.gp_para, xss(:,1:2), opt.regions, opt.level, 1, .8);
    opt.kdOBJ = KDTreeSearcher(xss(:,1:2));
%     opt.u = opt.method_obj.utility();
% % [idxInfo, ~] = knnsearch(kdOBJ,xstemp(1:opt.dim,:)');
% % f =-sum(opt.info(idxInfo,opt.dim+1));
end

ymu=zeros(size(xss(:,1)));
ys2=zeros(size(xss(:,1)));

if draw
    %--------------------------------------------------------%
    %%%%%%%%%%%%%%%%initialize figures%%%%%%%%%%%%%%%%%%%%%%%
    figure(1);set(gcf,'color','w');
    set(gcf, 'Position', [100, 400, 400, 400]);hold on;
    % h1=color_line3(xss(:,1), xss(:,2), xss(:,3),ymu,'.');
    h1= surface(X,Y,Z,reshape(ymu,size(X)));hold on
    shading interp
    opt.ceFig_optimal=  draw_path(traj(opt.z, opt), 'r', 2, 5);
    h1GP=scatter3(opt.gp.posGP(:,1),opt.gp.posGP(:,2),opt.gp.posGP(:,3),20,'filled','mo');
    colorbar
    axis equal tight
    axis([ DomainBounds.xmin DomainBounds.xmax DomainBounds.ymin DomainBounds.ymax])
    view(0,90)
    title('Predicted stiffness map and optimal trajectory')
    figure(2);set(gcf,'color','w');
    set(gcf, 'Position', [500, 400, 400, 400]);hold on;
    % h2=color_line3(xss(:,1), xss(:,2), xss(:,3),opt.gp.AF(:),'.') ;
    h2= surface(X,Y,Z,reshape(opt.gp.AF,size(X)));hold on
    shading interp
    axis equal tight
    opt.ceFig_candidate= draw_path(rand(3,3), 'b',3, 5);    %plot trajectories inside cem
    h2GP=scatter3(opt.gp.posGP(:,1),opt.gp.posGP(:,2),opt.gp.posGP(:,3),20,'filled','mo');
    axis equal
    axis([ DomainBounds.xmin DomainBounds.xmax DomainBounds.ymin DomainBounds.ymax])
    view(0,90)
    title('Acquisition function and candidate trajectories')
    figure(3);set(gcf,'color','w');
    set(gcf, 'Position', [900, 400, 400, 400]);hold on;
    % color_line3(xss(:,1), xss(:,2), xss(:,3),sGT,'.');
    h3= surface(X,Y,Z,reshape(sGT,size(X)));hold on
    shading interp
    axis equal tight
    axis([ DomainBounds.xmin DomainBounds.xmax DomainBounds.ymin DomainBounds.ymax])
    view(0,90)
    opt.ceFlag=0;


    % figure(4);set(gcf,'color','w');
    % set(gcf, 'Position', [1500, 400, 400, 400]);hold on;
    % h3 = color_line3(xss(:,1), xss(:,2), xss(:,3),ymu,'.');
    % title('Predicted stiffness map')
    % axis equal
    % axis([ DomainBounds.xmin DomainBounds.xmax DomainBounds.ymin DomainBounds.ymax])
end
%--------------------------------------------------------%

%% Receding-horizon trajectory planning

for k=1:queryLen %number of iterations
    opt.currentStage = k
    query_count=k-1;
%     if k==15
%         opt.planner = 1
%         display('Switch Planner ---> Utility Maximization')
%     end
%     if k == 1
%         opt.tf = 20;
%     else
%         opt.tf = 50;
%     end
    [opt,xs_traj] = gen_traj_CE(opt);
    
    %Execute the trajectory
    if(opt.planner==1)
        %Utility Maximization: pick the best sample along the trajectory
        [traj_utility] = EvaluateTrajUtility(xs_traj,opt);
        [~,idx_utility]=max(traj_utility);
        posGP_new=xs_traj(1:opt.dim,idx_utility)';
        traj_save=[traj_save,xs_traj(:,1:idx_utility)];
        %update the initial condition for trajectory
        xf=xs_traj(:,idx_utility);
    end
    
    if(opt.planner==2)
        %Ergodic coverage: pick multiple samples with fixed-sampling rate
        posGP_new=xs_traj(1:opt.dim,2:5:end)';
        traj_save=[traj_save,xs_traj];
        %update the initial condition for trajectory
        xf=xs_traj(:,end);
        opt = accumulate_CK(opt, xs_traj);% updates CK of the whole trajectory (for efficient calculation)
    end
    
    opt.xi = xf;
    %GPR update

    opt.gp.posGP=[opt.gp.posGP;posGP_new];
    kGP_new=EvaluateStiffnessKnn(posGP_new,xss,sGT);
    opt.gp.kGP= [opt.gp.kGP; kGP_new];
    %remove points that are too close to each other from the training set
    %[posGPfilt,kGPfilt]=consolidator(posGP,kGP,'max',2);%noise
    
    [ymu, ys2, fmu, fs2]= gp(opt.gp_para,opt.gp_model.inf, opt.gp_model.mean, opt.gp_model.cov, opt.gp_model.lik, opt.gp.posGP, opt.gp.kGP, xss);
    if opt.method==3
        opt.method_obj.update(posGP_new, kGP_new);
        lse.update(posGP_new(:, 1:2), kGP_new);
        region_measure_lse.update(posGP_new(:, 1:2), kGP_new);
        recall_lse(run_id, query_count+1) = (0+region_measure_lse.cumfound>0)'*region_outcome_gnd  / sum(region_outcome_gnd);
    end
    if opt.method==4
        opt.method_obj.update(posGP_new(:, 1:2), kGP_new);
        aas.update(posGP_new(:, 1:2), kGP_new);

        [~, ~, f, fs2] = aas.predict_points(x_gnd);
        recall_aas(run_id, query_count+1) = (0+aas.cumfound>0)'*region_outcome_gnd / sum(region_outcome_gnd);
    end
    
    if(opt.method==1)%Expected Improvement
        yEI=max(opt.gp.kGP);%current max
        opt.gp.AF = EI(yEI,ymu,ys2);
        ei.update(posGP_new(:, 1:2), kGP_new);
        region_measure_ei.update(posGP_new(:, 1:2), kGP_new);
        recall_ei(run_id, query_count+1) = (0+region_measure_ei.cumfound>0)'*region_outcome_gnd / sum(region_outcome_gnd);

    end
    if(opt.method==2)%Variance Reduction
        opt.gp.AF = ys2;
        unc.update(posGP_new(:, 1:2), kGP_new);
        region_measure_unc.update(posGP_new(:, 1:2), kGP_new);
        recall_unc(run_id, query_count+1) = (0+region_measure_unc.cumfound>0)'*region_outcome_gnd / sum(region_outcome_gnd);

    end
    if opt.method>=3
        opt.gp.AF = opt.method_obj.utility();
    end
    
    if(opt.planner==2)%Ergodic coverage
        opt.erg.mu=opt.gp.AF./sum(sum(opt.gp.AF));
        % Number of wave-numbers to be used
%         opt.erg.mu=flipud(opt.erg.mu);
%         opt.erg.mu=imrotate(opt.erg.mu,-90);
        [opt.erg.muk] = GetFourierCoeff(opt,X,Y);
        
    end
    if (draw)    
        set(h1,'XData', X,'YData', Y,'ZData',Z ,'CData', reshape(ymu,size(X)))
        set(h2,'XData', X,'YData',Y,'ZData',Z ,'CData', reshape(opt.gp.AF,size(X)))
        set(opt.ceFig_optimal,'XData', traj_save(1,:),'YData', traj_save(2,:),'ZData', traj_save(3,:));
        set(h1GP,'XData',opt.gp.posGP(:,1) ,'YData',opt.gp.posGP(:,2) ,'ZData', opt.gp.posGP(:,3));  %sampled points
        set(h2GP,'XData',opt.gp.posGP(:,1) ,'YData',opt.gp.posGP(:,2) ,'ZData', opt.gp.posGP(:,3));  %sampled points

        drawnow
    end
end
end
end
%%
zeros(size(mean(recall_aas, 1)))
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





 