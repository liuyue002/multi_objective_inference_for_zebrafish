%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load TDA scores, perform AABC, and visualize posterior 
%results for the iridophore transition parameter study.
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
% Note: the full dataset is too large to be hosted on Github. It can be
% obtained directly from the authors. The commented out code demonstrate
% how the posterior was obtained.

% load('../../example_results/inference_results/iridophore_results.mat','params_true','criteria','min_dists','num_criteria','num_nearest','num_presample','num_sample','num_trial','param_sample','presample_params','scaling','tda_dists','param_lbs','param_ubs','param_names');
% weights = 1./min_dists';
%% Combination of 192 objectives (recreates Fig.9)
% weighted_tda_dists = sum(tda_dists.*weights,2);
% num_accept=250;
% accept_eps=mink(weighted_tda_dists,num_accept+1);
% accept_eps=accept_eps(end);
% accept = weighted_tda_dists < accept_eps;
% accepted_params=param_sample(accept,:);

%% load the pre-computed posterior
load('../../example_results/inference_results/iridophore_results_slim.mat');
%% Plot
figure(Position=[10,10,1500,400]);
tiledlayout(1,6);
for param=1:6
    nexttile;
    histogram(accepted_params(:,param),param_lbs(param)-0.5:1:param_ubs(param)+0.5);
    hold on;
    xline(mean(accepted_params(:,param)),'r--',LineWidth=4);
    xline(params_true(param),'g--',LineWidth=4);
    xlabel(param_names{param},Interpreter="latex");
end
tl.TileSpacing="compact";
