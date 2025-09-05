% Inference for melanophore birth parameters
addpath('../tda');
%% Step 1: generate presamples, simulate the model, and TDA

savefolder='/path/to/repo/zebrafish_inference_public/example_results/melanophore/';
alpha_true = 3;
beta_true=3.5;
num_trial = 5;
tda_params=tda_default_params();
pimg_numpts = tda_params.pimg_numpts;
pland_kmax = tda_params.pland_kmax;
pland_numpts = tda_params.pland_numpts;
T=45;

% Ground truth simulations
% no nacre, because nacre lacks melanophores
for mutant = ["wt","pfeffer","shady"]
    for randomseed = 10001:10005
        melanophore_sim(alpha,beta,mutant,randomseed,savefolder);
    end
end

% Presample simulations. We use a regular lattice for the presamples
% instead of generating them randomly.
num_alpha = 21;
num_beta = 21;
num_presample = num_alpha * num_beta;
betas=linspace(0,10,num_beta);
alphas = linspace(0,10,num_alpha);
[Alphas,Betas] = meshgrid(alphas,betas);
presample_params = [Alphas(:),Betas(:)];

for mutant = ["wt","pfeffer","shady"] % no nacre
    for i=1:num_presample
        alpha = presample_params(i,1);
        beta = presample_params(i,2);
        for randomseed = 2024001:2024005
            % Note: this do all the simulations sequentially on the CPU,
            % which will take a very long time. You might want to write a
            % script to submit each simulation to a cluster to run them in
            % parallel, e.g. by wrapping melanophore_sim(...) in a slurm
            % wrapper, example given in melanophore_slurm.sh;
            % the implementation detail depends on your computing cluster
            melanophore_sim(alpha,beta,mutant,randomseed,savefolder);
        end
    end
end

%% Step 2: Gather TDA output

for mutant = ["wt","pfeffer","shady"]
    for celltype = ["M", "Xd", "Xl"]
        for dim=[0,1]
            for t = 1:45
                try
                    plands_true = cell(1,num_trial);
                    pland_true_avg=zeros(pland_kmax,pland_numpts);
                    for trial=1:num_trial
                        randomseed=10000+trial;
                        simname=sprintf('zebrafish_melanophores_%s_alpha=%.2f,beta=%.2f,rng=%d',mutant,alpha_true,beta_true,randomseed);
                        file=sprintf('%s/%s.mat',savefolder,simname);
                        data=load(file,'plands');
                        plands_true{trial} = data.plands.(celltype).(sprintf('dim%d',dim)){t};
                        if numel(plands_true{trial})==0
                            plands_true{trial} = zeros(pland_kmax,pland_numpts);
                        end
            			if size(plands_true{trial},1) > pland_kmax
            			    plands_true{trial} = plands_true{trial}(1:pland_kmax,:);
            			end
                        fprintf('%s : load OK\n', simname);
                        pland_true_avg = pland_true_avg + plands_true{trial};
                    end
                    pland_true_avg = pland_true_avg/num_trial;

                    plands = cell(num_presamples,num_trial);
                    avg_distances_pland = zeros(num_presamples,1);
                    for i=1:num_presamples
                        alpha = presample_params(i,1);
                        beta = presample_params(i,2);
                        pland_avg = zeros(pland_kmax,pland_numpts);
                        for trial=1:num_trial
                            randomseed=2024000+trial;
                            simname=sprintf('zebrafish_melanophores_%s_alpha=%.2f,beta=%.2f,rng=%d',mutant,alpha,beta,randomseed);
                            file=sprintf('%s/%s.mat',savefolder,simname);
                            data=load(file,'plands');
                            plands{i,trial} = data.plands.(celltype).(sprintf('dim%d',dim)){t};
                            if numel(plands{i,trial})==0
                                plands{i,trial} = zeros(pland_kmax,pland_numpts);
                            end
            			    if size(plands{i,trial},1) > pland_kmax
               			        plands{i,trial} = plands{i,trial}(1:pland_kmax,:);
            			    end
                            fprintf('%s : load OK\n', simname);
                            pland_avg = pland_avg + plands{i,trial};
                        end
                        pland_avg = pland_avg/num_trial;
                        avg_distances_pland(i) = sum((pland_avg-pland_true_avg).^2,'all')/(pland_kmax*pland_numpts);
                    end

                    avg_distances_pland_2d = reshape(avg_distances_pland,[num_beta,num_alpha]);
                    clear('data');

                    %% save
                    data=[]; % avoid saving it
                    save(sprintf('%s/melanophore_params_tda_%s_t=%02d_%s_dim%d.mat',savefolder,mutant,t,celltype,dim),'-v7.3');
                    clear('plands','avg_distances_pland','plands_true','pland_true_avg','avg_distances_pland_2d');

                catch exception
                    disp(exception);
                    fprintf('Melanophore_%s_t=%02d_%s_dim%d : ERROR, no output\n',mutant,t,celltype,dim);
                end
            end
        end
    end
end

%% Step 3: Decide which score function to use, and load relevant TDA output
tda_params=tda_default_params;
pland_kmax=tda_params.pland_kmax;
pland_numpts = tda_params.pland_numpts;

% uses 4 hand-picked score function
num_criteria = 4;
criteria = cell(num_criteria,1);
criteria{1} = struct('mutant','wt','t',26,'celltype','Xd','dim',1);
criteria{2} = struct('mutant','wt','t',45,'celltype','M','dim',1);
criteria{3} = struct('mutant','pfeffer','t',25,'celltype','M','dim',1);
criteria{4} = struct('mutant','shady','t',45,'celltype','Xd','dim',1);

% Alternatively, use a broad range of score functions
% num_criteria = 0;
% criteria = {};
% for i=24:45
%     num_criteria = num_criteria+1;
%     criteria{num_criteria} = struct('mutant','wt','t',i,'celltype','Xd','dim',1);
% end
% for i=24:45
%     num_criteria = num_criteria+1;
%     criteria{num_criteria} = struct('mutant','wt','t',i,'celltype','M','dim',1);
% end
% for i=20:45
%     num_criteria = num_criteria+1;
%     criteria{num_criteria} = struct('mutant','pfeffer','t',i,'celltype','M','dim',1);
% end
% for i=36:45
%     num_criteria = num_criteria+1;
%     criteria{num_criteria} = struct('mutant','shady','t',i,'celltype','Xd','dim',1);
% end

plands_presample = cell(num_criteria,1);
plands_true_avg = cell(num_criteria,1);
plands_true_differences = zeros(num_criteria,num_trial);
min_dists = zeros(num_criteria,1);

for i=1:num_criteria
    file=sprintf('%s/melanophore_params_tda_%s_t=%02d_%s_dim%d.mat',savefolder,criteria{i}.mutant,criteria{i}.t,criteria{i}.celltype,criteria{i}.dim);
    data=load(file,'pland_true_avg','plands','avg_distances_pland_2d','plands_true');
    plands_presample{i} = data.plands;
    plands_true_avg{i} = data.pland_true_avg(1:pland_kmax,:);
    min_dists(i) = min(data.avg_distances_pland_2d,[],'all');
    for trial=1:num_trial
        plands_true_differences(i,trial) = sum((data.plands_true{trial}-plands_true_avg{i}).^2,'all')/(pland_kmax*pland_numpts);
    end
end

%% Step 4: AABC: Compute score functions

num_sample=100000;
num_nearest = 5;
param_sample=rand(num_sample,2)*10; % prior ~ unif(0,10)
tda_dists=zeros(num_sample,num_criteria);

tic;
for i=1:num_sample
    sample=param_sample(i,:);
    param_dist = sqrt(sum((presample_params./scaling-sample./scaling).^2,2));
    dist_k1 = mink(param_dist,num_nearest+1);
    dist_k1 = dist_k1(end);
    idx_k = find(param_dist==dist_k1,1);
    ws = (3/4)*(1/dist_k1)*(1- (param_dist/dist_k1).^2 ).*(param_dist<dist_k1); % weights for neighbors
    phi = drchrnd(ws,1);
    sample_idx = randsample(num_presamples,num_trial,true,phi);
    for j=1:num_criteria
        avg_pland=zeros(pland_kmax,pland_numpts);
        for k=1:num_trial
            idx = randi(num_trial);
            avg_pland = avg_pland + plands_presample{j}{sample_idx(k),idx}(1:pland_kmax,:);
        end
        avg_pland = avg_pland/num_trial;
        tda_dists(i,j) = sum((avg_pland-plands_true_avg{j}).^2,'all')/(pland_kmax*pland_numpts);
    end
end
fprintf('Computing TDA distance for %d samples took %f sec\n',num_sample,toc);

%% Step 5: AABC: Acceptance by weighted sum
% weight each criteria by the min dist achived by a presample point
weights = 1./min_dists';
weighted_tda_dists = sum(tda_dists.*weights,2);
accept_proportion=0.01; % important hyper-parameter to tune
accept_eps=mink(weighted_tda_dists,round(accept_proportion*num_sample)+1);
accept_eps=accept_eps(end,:);

accept = weighted_tda_dists < accept_eps;
num_accept = sum(accept);
accepted_params=param_sample(accept,:);

%% Step 6: Visualise posterior distribution and save

fig_posterior=figure;
hold on
plot(accepted_params(:,1),accepted_params(:,2),'.');
plot(mean(accepted_params(:,1)),mean(accepted_params(:,2)),'.r',MarkerSize=20);
plot(3,3.5,'.g',MarkerSize=20);
xlim([0,10]);
ylim([0,10]);
xlabel('$\alpha$',interpreter="latex");
ylabel('$\beta$',interpreter="latex");
title('AABC posterior');

saveas(fig_posterior,sprintf('%s/melanophore_posterior.png',savefolder));
saveas(fig_posterior,sprintf('%s/melanophore_posterior.eps',savefolder),'epsc');
save(sprintf('%s/melanophore_aabc_results',savefolder),'-v7.3');