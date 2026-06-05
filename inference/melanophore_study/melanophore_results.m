%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load TDA scores, perform AABC, and visualize posterior 
%results for the melanophore differentiation parameter study.
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

%% load data
load('../../example_results/inference_results/melanophore_inference_results.mat','Alphas','Betas','alphas','betas','criteria','min_dists','num_criteria','num_nearest','num_presamples','num_sample','num_trial','param_sample','presample_params','scaling','tda_dists');
weights = 1./min_dists';

%% Single objective (recreates Fig.6)

% score_indx select the score function corresponding to each subfigure
for score_indx = [44,3,80,50,23]
    score_fn=criteria{score_indx};
    accept_eps=1.2;
    accept = tda_dists(:,score_indx) < accept_eps;
    num_accept = sum(accept);
    accepted_params=param_sample(accept,:);

    figure;
    hold on
    plot(accepted_params(:,1),accepted_params(:,2),'.');
    plot(mean(accepted_params(:,1)),mean(accepted_params(:,2)),'.r',MarkerSize=20);
    plot(3,3.5,'.g',MarkerSize=20);
    xlim([0,10]);
    ylim([0,10]);
    title(sprintf('%s, %s, t=%d, dim=1, mean=$(%.2f, %.2f)$, accepted %d/%d',score_fn.mutant,score_fn.celltype,score_fn.t,mean(accepted_params),num_accept,num_sample),Interpreter="latex");
    xlabel('$\alpha$',Interpreter='latex');
    ylabel('$\beta$',Interpreter='latex');
end

%% Combination of 4 objectives (recreates Fig.7 bottom left)

score_indx = [44,3,80,50];
tda_dists4=tda_dists(:,score_indx);
weighted_tda_dists4 = sum(tda_dists4.*weights(score_indx),2);
accept_proportion=0.01;
accept_eps=mink(weighted_tda_dists4,round(accept_proportion*num_sample)+1);
accept_eps=accept_eps(end,:);
accept = weighted_tda_dists4 < accept_eps;
num_accept = sum(accept);
accepted_params=param_sample(accept,:);

figure;
hold on
plot(accepted_params(:,1),accepted_params(:,2),'.');
plot(mean(accepted_params(:,1)),mean(accepted_params(:,2)),'.r',MarkerSize=20);
plot(3,3.5,'.g',MarkerSize=20);
xlim([0,10]);
ylim([0,10]);
xlabel('$\alpha$',interpreter="latex");
ylabel('$\beta$',interpreter="latex");
title(sprintf('4 objective weighted sum, mean=$(%.2f,%.2f)$',mean(accepted_params)),Interpreter="latex");
xticks([0,10]);
yticks([0,10]);

% note the resulting figure may be slightly different from the paper due to
% different random number generation used when sampling neighbours

%% Combination of 80 objectives (recreates Fig.7 bottom right)
weighted_tda_dists = sum(tda_dists.*weights,2);
accept_proportion=0.01;
accept_eps=mink(weighted_tda_dists,round(accept_proportion*num_sample)+1);
accept_eps=accept_eps(end,:);
accept = weighted_tda_dists < accept_eps;
num_accept = sum(accept);
accepted_params=param_sample(accept,:);

figure;
hold on
plot(accepted_params(:,1),accepted_params(:,2),'.');
plot(mean(accepted_params(:,1)),mean(accepted_params(:,2)),'.r',MarkerSize=20);
plot(3,3.5,'.g',MarkerSize=20);
xlim([0,10]);
ylim([0,10]);
xlabel('$\alpha$',interpreter="latex");
ylabel('$\beta$',interpreter="latex");
title(sprintf('80 objective weighted sum, mean=$(%.2f,%.2f)$',mean(accepted_params)),Interpreter="latex");
xticks([0,10]);
yticks([0,10]);