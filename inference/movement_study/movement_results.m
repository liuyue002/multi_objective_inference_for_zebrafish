%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load TDA scores, perform AABC, and visualize posterior 
%results for the movement parameters study.
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
load('D:\zebrafish_storage\maclaptop\movement_20241111/post-analysis/aabc_29crit.mat','RIDIDS','RILILS','Ridids','Rilils','Ridid_true','Rilil_true','criteria','min_dists','num_criteria','num_nearest','num_presample','num_sample','num_trial','param_sample','presample_params','scaling','tda_dists');
weights = 1./min_dists';
%% Single objective (recreates Fig.8 a-d)

% score_indx select the score function corresponding to each subfigure
for score_indx = [11,6,13,20]
    score_fn=criteria{score_indx};
    accept_eps=0.8;
    accept = tda_dists(:,score_indx) < accept_eps;
    num_accept = sum(accept);
    accepted_params=param_sample(accept,:);

    figure;
    hold on
    plot(accepted_params(:,1),accepted_params(:,2),'.');
    plot(mean(accepted_params(:,1)),mean(accepted_params(:,2)),'.r',MarkerSize=20);
    plot(Ridid_true,Rilil_true,'.g',MarkerSize=20);
    xlim([0,20]);
    ylim([0,80]);
    title(sprintf('%s, %s, t=%d, dim=1, mean=$(%.2f, %.2f)$, accepted %d/%d',score_fn.mutant,score_fn.celltype,score_fn.t,mean(accepted_params),num_accept,num_sample),Interpreter="latex");
    xlabel('$R_{I^d,I^d}$',Interpreter='latex');
    ylabel('$R_{I^l,I^l}$',Interpreter='latex');
end

%% Combination of 29 objectives (recreates Fig.8e)
weighted_tda_dists = sum(tda_dists.*weights,2);
accept_eps=num_criteria*7;
accept = weighted_tda_dists < accept_eps;
num_accept = sum(accept);
accepted_params=param_sample(accept,:);

figure;
hold on
plot(accepted_params(:,1),accepted_params(:,2),'.');
plot(mean(accepted_params(:,1)),mean(accepted_params(:,2)),'.r',MarkerSize=20);
plot(Ridid_true,Rilil_true,'.g',MarkerSize=20);
xlim([0,20]);
ylim([0,80]);
xlabel('$\alpha$',interpreter="latex");
ylabel('$\beta$',interpreter="latex");
title(sprintf('29 objective weighted sum, mean=$(%.2f,%.2f)$',mean(accepted_params)),Interpreter="latex");
xlabel('$R_{I^d,I^d}$',Interpreter='latex');
ylabel('$R_{I^l,I^l}$',Interpreter='latex');