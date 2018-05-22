%Authors: Hadi Salman, Elif Ayvali, and Yifei Ma
%% Initialization
clc;close all;clear;
addpath(genpath(pwd))
colormap hot
name_exp = ['exp_1___', date];
mkdir(['results/' name_exp]);
%write something about this experiment
fid = fopen( ['results/' name_exp '/readme.txt'], 'wt' );
readme='number of Stages = 60....probes per stage = 2 .... level = 1....ell = 6... beta_t = .5...eps_band=0.1 ';
fprintf(fid,readme);
fclose(fid);
mkdir(['results/' name_exp '/ground_truth'])
mkdir(['results/' name_exp '/aas/estimated'])
mkdir(['results/' name_exp '/aas/acquisition'])
mkdir(['results/' name_exp '/lse/estimated'])
mkdir(['results/' name_exp '/lse/acquisition'])
mkdir(['results/' name_exp '/unc/estimated'])
mkdir(['results/' name_exp '/unc/acquisition'])
mkdir(['results/' name_exp '/ei/estimated'])
mkdir(['results/' name_exp '/ei/acquisition'])
mkdir(['results/' name_exp '/rnd/estimated'])
mkdir(['results/' name_exp '/rnd/acquisition'])

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

%%%%%%%%%%%%%%%%Import Data%%%%%%%%%%%%%%%%%%%%%
% load('1tumor_1artery.mat');
load('2tumors_experiment.mat');
xss=x_gnd;%3D grid pos
x1 = reshape(xss(:,1),51,51);
x2 = reshape(xss(:,2),51,51);
x3 = reshape(xss(:,3),51,51);

sGT=y_gnd;%Stiffness ground truth
% sGT=y_gnd/max(max(y_gnd));%Stiffness ground truth
%(3d point, scalar stiffness at the point)
datafull=[xss,sGT(:)];
%--------------------------------------------------------%
%%%%%%%%%%%%%%%%Setting domain bounds%%%%%%%%%%%%%%%%%%%%%
DomainBounds.xmin = min(xss(:,1));
DomainBounds.xmax = max(xss(:,1));
DomainBounds.ymin = min(xss(:,2));
DomainBounds.ymax = max(xss(:,2));
Lx = DomainBounds.xmax - DomainBounds.xmin;
Ly = DomainBounds.ymax - DomainBounds.ymin;
%--------------------------------------------------------%
for alg=1:4
%%%%%%%%%%Initialize Prior Acquisition Function%%%%%%%%%%%%%
%Needs to be defined on a 2D grid
opt = [];%initialize options
opt.gp.AF=zeros(sqrt(length(xss(:,1))),sqrt(length(xss(:,2)))); 
%--------------------------------------------------------%
%%%%%%%%%%%%%%initialize robot position%%%%%%%%%%%%%%%%%%%%
%Get te index of the position where AF is max
% [~,idxinit]=max(opt.gp.AF(:));
idxinit=800;
%--------------------------------------------------------%
%%%%%%%%%%%%%%%%initialize CE planner%%%%%%%%%%%%%%%%%%%%%%%
[opt] = initialize_gen_traj_CE(opt,DomainBounds);
opt.stages = 30;
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

if(opt.method==3)
    display('Method: Level Set Estimation')
end

if(opt.method==4)
    display('Method: Active Area Search')
end

opt.xi = [xss(idxinit,:)'; 0.1*pi];%[pos;orientation]
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
    opt.erg.mu=flipud(opt.erg.mu);
    opt.erg.mu=imrotate(opt.erg.mu,-90);
    [opt.erg.muk] = GetFourierCoeff(opt,x1,x2);
end
%--------------------------------------------------------%
%%%%%%%%%%%%%%%%initialize GP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GPRsave.ymusave=[];GPRsave.ys2save=[]; GPRsave.AFsave=[];
opt.gp.posGP=[];opt.gp.kGP=[];sample_rate=5;
%GP parameters
opt.gp_model = struct('inf',@infExact, 'mean', @meanZero, 'cov', @covSEiso, 'lik', @likGauss);
sn = 0.01; ell =6; sf = sqrt(1);Ncg=30;
opt.gp_para.lik = log(sn); opt.gp_para.cov = log([ell; sf]);
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
    [region_grid_x, region_grid_y] = meshgrid(0:5:75, 0:5:80);
    opt.regions = [region_grid_x(:), region_grid_x(:)+5, region_grid_y(:), region_grid_y(:)+5];
%     res = 20;
%     [region_grid_x, region_grid_y] = meshgrid(linspace(DomainBounds.xmin,DomainBounds.xmax,res), linspace(DomainBounds.ymin,DomainBounds.ymax,res));
%     region_grid_x = region_grid_x(:,1:end-1);
%     region_grid_y = region_grid_y(1:end-1,:);
%     opt.regions = [region_grid_x(:), region_grid_x(:)+Lx/res, region_grid_y(:), region_grid_y(:)+Ly/res];
    opt.method_obj = ActiveAreaSearch(opt.gp_model, opt.gp_para, xss(:,1:2), opt.regions, opt.level, 1, .8);
    opt.kdOBJ = KDTreeSearcher(xss(:,1:2));
%     opt.u = opt.method_obj.utility();
% % [idxInfo, ~] = knnsearch(kdOBJ,xstemp(1:opt.dim,:)');
% % f =-sum(opt.info(idxInfo,opt.dim+1));
end

ymu=zeros(size(xss(:,1)));
ys2=zeros(size(xss(:,1)));
GPRsave.ymusave(1,:)=ymu;
GPRsave.ys2save(1,:)=ys2;
GPRsave.AFsave(1,:)=opt.gp.AF(:);


%--------------------------------------------------------%
%%%%%%%%%%%%%%%%initialize figures%%%%%%%%%%%%%%%%%%%%%%%
figure(1);set(gcf,'color','w');clf
set(gcf, 'Position', [100, 400, 400, 400]);hold on;
% h1=color_line3(xss(:,1), xss(:,2), xss(:,3),ymu,'.');
h1= surface(x1,x2,x3,reshape(ymu,size(x1)));hold on
shading interp
opt.ceFig_optimal=  draw_path(traj(opt.z, opt), 'c', .2, 1);
h1GP=scatter3(opt.gp.posGP(:,1),opt.gp.posGP(:,2),opt.gp.posGP(:,3),10,'filled','k');
colorbar
axis equal tight
axis([ DomainBounds.xmin DomainBounds.xmax DomainBounds.ymin DomainBounds.ymax])
view(0,90)
title('Predicted stiffness map and optimal trajectory')
figure(2);set(gcf,'color','w');clf
set(gcf, 'Position', [500, 400, 400, 400]);hold on;
% h2=color_line3(xss(:,1), xss(:,2), xss(:,3),opt.gp.AF(:),'.') ;
h2= surface(x1,x2,x3,reshape(opt.gp.AF,size(x1)));hold on
shading interp
axis equal tight
opt.ceFig_candidate= draw_path(rand(3,3), 'b',3, 5);    %plot trajectories inside cem
h2GP=scatter3(opt.gp.posGP(:,1),opt.gp.posGP(:,2),opt.gp.posGP(:,3),20,'filled','mo');
axis equal
axis([ DomainBounds.xmin DomainBounds.xmax DomainBounds.ymin DomainBounds.ymax])
view(0,90)
title('Acquisition function and candidate trajectories')
figure(3);set(gcf,'color','w');clf
set(gcf, 'Position', [900, 400, 400, 400]);hold on;
% color_line3(xss(:,1), xss(:,2), xss(:,3),sGT,'.');
h3= surface(x1,x2,x3,reshape(sGT,size(x1)));hold on
shading interp
colormap hot
axis equal tight
axis([ DomainBounds.xmin DomainBounds.xmax DomainBounds.ymin DomainBounds.ymax])
view(0,90)
opt.ceFlag=0;
%--------------------------------------------------------%

%% Receding-horizon trajectory planning

for k=1:opt.stages %number of iterations
    opt.currentStage = k
    if k==1
       opt.tf = 5;
    else
        opt.tf = 20;
    end
    
    [opt,xs_traj] = gen_traj_CE(opt);
    %Execute the trajectory
    if(opt.planner==1)
        %Utility Maximization: pick the best sample along the trajectory
        [traj_utility] = EvaluateTrajUtility(xs_traj,opt);
        [~,idx_utility]=max(traj_utility);
%         posGP_new=xs_traj(1:opt.dim,idx_utility)';
        traj_save=[traj_save,xs_traj(:,1:idx_utility)];
        xs_traj = xs_traj(:,1:idx_utility)
        posGP_new=xs_traj(1:opt.dim,end:-10:1)';
        %update the initial condition for trajectory
        xf=xs_traj(:,end);
    end
    
    if(opt.planner==2)
        %Ergodic coverage: pick multiple samples with fixed-sampling rate
        posGP_new=xs_traj(1:opt.dim,1:10:end)';
        traj_save=[traj_save,xs_traj];
        %update the initial condition for trajectory
        xf=xs_traj(:,end);
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
        for i=1:length(kGP_new)
            opt.method_obj.update(posGP_new(i,:), kGP_new(i));
        end
    end
    if opt.method==4
        for i=1:length(kGP_new)
            opt.method_obj.update(posGP_new(i,1:2), kGP_new(i));
        end
    end
    
    if(opt.method==1)%Expected Improvement
        yEI=max(opt.gp.kGP);%current max
        opt.gp.AF = EI(yEI,ymu,ys2);
    end
    if(opt.method==2)%Variance Reduction
        opt.gp.AF = ys2;
    end
    if opt.method>=3
        opt.gp.AF = opt.method_obj.utility();
    end
    
    if(opt.planner==2)%Ergodic coverage
        opt.erg.mu=opt.gp.AF./sum(sum(opt.gp.AF));
        % Number of wave-numbers to be used
        opt.erg.mu=flipud(opt.erg.mu);
        opt.erg.mu=imrotate(opt.erg.mu,-90);
        [opt.erg.muk] = GetFourierCoeff(opt,x1,x2);
    end
    
    
    set(h1,'XData', x1,'YData', x2,'ZData',x3 ,'CData', reshape(ymu,size(x1)))
    set(h2,'XData', x1,'YData',x2,'ZData',x3 ,'CData', reshape(opt.gp.AF,size(x1)))
    set(opt.ceFig_optimal,'XData', traj_save(1,:),'YData', traj_save(2,:),'ZData', traj_save(3,:));
    set(h1GP,'XData',opt.gp.posGP(:,1) ,'YData',opt.gp.posGP(:,2) ,'ZData', opt.gp.posGP(:,3));  %sampled points
    set(h2GP,'XData',opt.gp.posGP(:,1) ,'YData',opt.gp.posGP(:,2) ,'ZData', opt.gp.posGP(:,3));  %sampled points
    drawnow
    
    figure(1)
    if(alg==1)
        saveas(gcf,['results/' name_exp '/ei/estimated/ei_',int2str(k)])
    elseif(alg==2)
        saveas(gcf,['results/' name_exp '/unc/estimated/unc_',int2str(k)])
    elseif(alg==3)
        saveas(gcf,['results/' name_exp '/lse/estimated/lse_',int2str(k)])
    elseif(alg==4)
        saveas(gcf,['results/' name_exp '/aas/estimated/aas_',int2str(k)])
    end    
    figure(2)
    if(alg==1)
        saveas(gcf,['results/' name_exp '/ei/acquisition/ei_',int2str(k)])
    elseif(alg==2)
        saveas(gcf,['results/' name_exp '/unc/acquisition/unc_',int2str(k)])
    elseif(alg==3)
        saveas(gcf,['results/' name_exp '/lse/acquisition/lse_',int2str(k)])
    elseif(alg==4)
        saveas(gcf,['results/' name_exp '/aas/acquisition/aas_',int2str(k)])
    end    

        
    
    GPRsave.ymusave(opt.stages+1,:)=ymu;
    GPRsave.ys2save(opt.stages+1,:)=ys2;
    GPRsave.AFsave(opt.stages+1,:)=opt.gp.AF(:);
    
end
end
save(['results/' name_exp '/data'])
