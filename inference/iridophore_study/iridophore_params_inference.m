%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inference pipeline for iridophore transition parameters.
%
%
%Copyright (C) Yue Liu, Alexandria Volkening, 2026
%
%This program is free software: you can redistribute it and/or 
%modify it under the terms of the GNU General Public License as 
%published by the Free Software Foundation, version 3.
%
%This program is distributed in the hope that it will be useful, 
%but WITHOUT ANY WARRANTY; without even the implied warranty of 
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
%See the GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public 
%License along with this program. If not, see 
%https://www.gnu.org/licenses/.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../tda');
%% Step 1: generate presamples, simulate the model, and TDA

savefolder='/path/to/repo/zebrafish_inference_public/example_results/iridophore/';
params_true = [3,9,3,3,5,2];
num_trial = 5;
tda_params=tda_default_params();
pimg_numpts = tda_params.pimg_numpts;
pland_kmax = tda_params.pland_kmax;
pland_numpts = tda_params.pland_numpts;
T=45;

% Ground truth simulations
% no shady, because shady lacks iridophores
for mutant = ["wt","pfeffer","nacre"]
    for randomseed = 10001:10005
        iridophore_sim(params_true,mutant,randomseed,savefolder);
    end
end

% Presample simulations. 
% uniform prior: c=[0-15], d=[0-30], e=[-1,15], f=[-1,15], g=[-1,30], h=[0-15]
num_presample=4750;
param_lbs=[0,0,-1,-1,-1,0];
param_ubs=[15,30,15,15,30,15];
presample_params=zeros(num_sample,6);
for i=1:6
    presample_params(:,i) = randi(param_ubs(i)-param_lbs(i)+1,num_sample,1)-1+param_lbs(i);
end

for mutant = ["wt","pfeffer","nacre"] % no shady
    for i=1:num_presample
		irid_params=presample_params(i,:);
        for randomseed = 2024001:2024005
            % Note: this do all the simulations sequentially on the CPU,
            % which will take a very long time. You might want to write a
            % script to submit each simulation to a cluster to run them in
            % parallel, e.g. by wrapping iridophore_sim(...) in a slurm
            % wrapper, see the melanophore study for example.
            iridophore_sim(irid_params,mutant,randomseed,savefolder);
        end
    end
end

%% Step 2: Gather TDA output

for mutant = ["wt","pfeffer","nacre"]
    for celltype = ["M", "Xd", "Xl"]
        for dim=[0,1]
            for t = 1:45
                try
                    plands_true = cell(1,num_trial);
                    pland_true_avg=zeros(pland_kmax,pland_numpts);
                    for trial=1:num_trial
                        randomseed=10000+trial;
                        simname=sprintf('zebrafish_iridophore_[%02d,%02d,%02d,%02d,%02d,%02d]_%s,rng=%d',params_true,mutant,randomseed);
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
                        pland_avg = zeros(pland_kmax,pland_numpts);
						irid_params = presample_params(i,:);
                        for trial=1:num_trial
                            randomseed=2024000+trial;
                            simname=sprintf('zebrafish_iridophore_[%02d,%02d,%02d,%02d,%02d,%02d]_%s,rng=%d',irid_params,mutant,randomseed);
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
                    clear('data');

                    %% save
                    data=[]; % avoid saving it
                    save(sprintf('%s/iridophore_params_tda_%s_t=%02d_%s_dim%d.mat',savefolder,mutant,t,celltype,dim),'-v7.3');
                    clear('plands','avg_distances_pland','plands_true','pland_true_avg');

                catch exception
                    disp(exception);
                    fprintf('Iridophore_%s_t=%02d_%s_dim%d : ERROR, no output\n',mutant,t,celltype,dim);
                end
            end
        end
    end
end

%% Step 3: Decide which score function to use, and load relevant TDA output
tda_params=tda_default_params;
pland_kmax=tda_params.pland_kmax;
pland_numpts = tda_params.pland_numpts;

% uses basically every sensible score function
num_criteria = 0;
criteria = {};
for i=24:45
    num_criteria = num_criteria+1;
    criteria{num_criteria} = struct('mutant','wt','t',i,'celltype','M','dim',0);
end
for i=24:45
    num_criteria = num_criteria+1;
    criteria{num_criteria} = struct('mutant','wt','t',i,'celltype','M','dim',1);
end
for i=24:45
    num_criteria = num_criteria+1;
    criteria{num_criteria} = struct('mutant','wt','t',i,'celltype','Xd','dim',0);
end
for i=24:45
    num_criteria = num_criteria+1;
    criteria{num_criteria} = struct('mutant','wt','t',i,'celltype','Xd','dim',1);
end
for i=20:45
    num_criteria = num_criteria+1;
    criteria{num_criteria} = struct('mutant','pfeffer','t',i,'celltype','M','dim',0);
end
for i=20:45
    num_criteria = num_criteria+1;
    criteria{num_criteria} = struct('mutant','pfeffer','t',i,'celltype','M','dim',1);
end
for i=20:45
    num_criteria = num_criteria+1;
    criteria{num_criteria} = struct('mutant','nacre','t',i,'celltype','Xd','dim',0);
end
for i=20:45
    num_criteria = num_criteria+1;
    criteria{num_criteria} = struct('mutant','nacre','t',i,'celltype','Xd','dim',1);
end

plands_presample = cell(num_criteria,1);
plands_true_avg = cell(num_criteria,1);
plands_true_differences = zeros(num_criteria,num_trial);
min_dists = zeros(num_criteria,1);

for i=1:num_criteria
    file=sprintf('%s/iridophore_params_tda_%s_t=%02d_%s_dim%d.mat',savefolder,criteria{i}.mutant,criteria{i}.t,criteria{i}.celltype,criteria{i}.dim);
    data=load(file,'pland_true_avg','plands','avg_distances_pland','plands_true');
    plands_presample{i} = data.plands;
    plands_true_avg{i} = data.pland_true_avg(1:pland_kmax,:);
    min_dists(i) = min(data.avg_distances_pland);
    for trial=1:num_trial
        plands_true_differences(i,trial) = sum((data.plands_true{trial}-plands_true_avg{i}).^2,'all')/(pland_kmax*pland_numpts);
    end
end

%% Step 4: AABC: Compute score functions
rng(12345);
num_sample=1000000;
num_nearest = 5;
param_sample=zeros(num_sample,6);
for i=1:6
    param_sample(:,i) = randi(param_ubs(i)-param_lbs(i)+1,num_sample,1)-1+param_lbs(i);
end
% add the presample params to samples
param_sample=[presample_params;param_sample];
num_sample=size(param_sample,1);
tda_dists=zeros(num_sample,num_criteria);
scaling=[1,1,1,1,1,1];

tic;
for i=1:num_sample
    sample=param_sample(i,:);
    param_dist = sqrt(sum((presample_params./scaling-sample./scaling).^2,2));
    dist_k1 = mink(param_dist,num_nearest+1);
    % deal with the case that the num_nearest+1 presamples all have same
    % dist (occasionally happens since params are integer-valued)
    if dist_k1(end) == dist_k1(end-1)
        sorted_dist = sort(param_dist);
        idx = find(sorted_dist > dist_k1(end),true,'first');
        num_nearest = idx-1;
        dist_k1 = mink(param_dist,num_nearest+1);
    end
	dist_k1 = dist_k1(end);
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
tl=tiledlayout(1,6);
for param=1:6
    nexttile;
    histogram(accepted_params(:,param),param_lbs(param)-0.5:1:param_ubs(param)+0.5);
    hold on;
    xline(mean(accepted_params(:,param)),'r--');
    xline(params_true(param),'g--');
    xlabel(param_names{param});
end
sgtitle('AABC posterior');

saveas(fig_posterior,sprintf('%s/iridophore_posterior.png',savefolder));
saveas(fig_posterior,sprintf('%s/iridophore_posterior.eps',savefolder),'epsc');
save(sprintf('%s/iridophore_aabc_results',savefolder),'-v7.3');