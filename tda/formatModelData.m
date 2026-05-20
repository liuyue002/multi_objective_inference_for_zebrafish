%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Utility function for computing pairwise distance and
%formatting output data from a simulation of the zebrafish
%model, in preparation for TDA computation.
%
%
%Copyright (C) Yue Liu, Alexandria Volkening, 2026
%Adapted from Cleveland, Zhu, Sandstede, & Volkening (2023)
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

function [output_file,formatted_data,pw_distances] = formatModelData(savefolder,simname,t,save_files,Ycrop)
%Format the output of the zebrafish simulation to prepare for TDA
%calculatios. This file is adapted from the code by Cleveland et al.
% savefolder: path to folder for saving results
% simname: name of the simulation for saving results
% t: time of interest 
% (unit: days, t=1 correspond to 21 days post fertilization, t=45 ~ 65 dpf)
% save_files: whether to save intermediate files
% crop: How much of the domain in y-direction is cropped out 
% (default in Liu & Volkening and Cleveland et al is 0.1)
addpath('../TDACodeOriginal');

%%% Specifying other formatting patterns
Ycutoff = Ycrop;

load(strcat(savefolder,simname),'cellsM', 'cellsXsn', 'cellsXc', 'cellsIl', 'cellsId','boundaryX','boundaryY');
end_dpf = 20+size(boundaryX,2);       % (in days post fertilization)
boundaryX = boundaryX(end);      % extracting the length of the domain (and height of the domain next)
boundaryY = boundaryY(end);      % simulation starts at 21 dpf, and time step 46 corresponds to 66 dpf
heightmm = boundaryY/1000;      % converting untis from um to mm
SL = 12.63;                     % in mm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loading data associated with the agent-based model by Volkening
%%% and Sandstede (2018). Extra space in the matrices is marked as
%%% -40000. 
cellsM = cellsM(cellsM(:,1,t+1) > -39000,:,t+1);        % melanophores
cellsXl = cellsXsn(cellsXsn(:,1,t+1) > -39000,:,t+1);   % loose xanthophores
cellsXd = cellsXc(cellsXc(:,1,t+1) > -39000,:,t+1);     % dense xanthophores
cellsIl = cellsIl(cellsIl(:,1,t+1) > -39000,:,t+1);     % loose iridophores
cellsId = cellsId(cellsId(:,1,t+1) > -39000,:,t+1);     % dense iridophores
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Gathering the number of cells of each type
numM = size(cellsM,1);
numXl = size(cellsXl,1);
numXd = size(cellsXd,1);
numIl = size(cellsIl,1);
numId = size(cellsId,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% For consistency with the cellular automaton model by Owen, 
%%% Kelsh, and Yates (2020), we assign a number label to each cell 
%%% type.
identifierM = ones(numM,1);         % identifier 1 marks melanophores
identifierXl = 2*ones(numXl,1);     % identifier 2 marks loose xanthophores
identifierXd = 4*ones(numXd,1);     % identifier 4 marks dense xanthophores
identifierIl = 5*ones(numIl,1);     % identifier 5 marks loose iridophores
identifierId = 6*ones(numId,1);     % identifier 6 marks dense iridophores

%%% Combining cell positions and their identifiers
fullM = [cellsM, identifierM];
fullXsn = [cellsXl, identifierXl];
fullXc = [cellsXd, identifierXd];
fullIl = [cellsIl, identifierIl];
fullId = [cellsId, identifierId];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Defining one cells_array variable with data for all cell types
cells_array = [fullM; fullXsn; fullXc; fullIl; fullId];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Cropping patterns and gathering cropped cells. We remove cells from
%%% the top Ycutoff% and bottom Ycutoff% of the domain
[cells_array_restricted] = cell_domain_restriction(cells_array,boundaryY,Ycutoff);

cells_mel_u = cells_array_restricted(cells_array_restricted(:,3)==1,:);   % unclean, cropped melanophores
cells_iriL = cells_array_restricted(cells_array_restricted(:,3)==5,:);    % unclean, cropped loose iridophores
cells_iriD = cells_array_restricted(cells_array_restricted(:,3)==6,:);    % unclean, cropped dense iridophores
cells_xanD_u = cells_array_restricted(cells_array_restricted(:,3)==4,:);  % unclean, cropped dense xanthophores
cells_xanL_u = cells_array_restricted(cells_array_restricted(:,3)==2,:);  % unclean, cropped loose xanthophores

%%% Updating the cells_array variable with our formatted data for all 
%%% types
cells_array = [cells_mel_u; cells_xanL_u ; cells_xanD_u; cells_iriL; cells_iriD];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% saving
if save_files
    output_file = sprintf('%s/%s_formatted_t=%d.mat',savefolder,simname,t);
    save(output_file,'cells_array', 'end_dpf','boundaryX','boundaryY', 'heightmm','SL');
else
    output_file='';
end
formatted_data.cells_array = cells_array;
formatted_data.end_dpf = end_dpf;
formatted_data.boundaryX = boundaryX;
formatted_data.boundaryY = boundaryY;
formatted_data.heightmm = heightmm;
formatted_data.SL = SL;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clean data step

preBCnum_nbrs = 10;         % number of nearest neighbors (excluding self) as above
preBCthresh_mel = 220;      % as above (in um)
preBCthresh_xanL = 170;     % in um
preBCthresh_xanD = 130;     % in um

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Gathering cell positions (that have already been croppped and
%%% formatted in Step 1)
cells_mel_u = cells_array(cells_array(:,3)==1,1:2);   % unclean, cropped melanophores
cells_iriL = cells_array(cells_array(:,3)==5,1:2);    % unclean, cropped loose iridophores
cells_iriD = cells_array(cells_array(:,3)==6,1:2);    % unclean, cropped dense iridophores
cells_xanD_u = cells_array(cells_array(:,3)==4,1:2);  % unclean, cropped dense xanthophores
cells_xanL_u = cells_array(cells_array(:,3)==2,1:2);  % unclean, cropped loose xanthophores
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Cleaning data to remove outlier cells (and tracking the number of
%%% cells cleaned pre computing barcodes)
[logical_clean_mel, logical_clean_xanD, logical_clean_xanL, ~, ~, ~] = ...
    cell_domain_cleaning(cells_mel_u, cells_xanD_u, cells_xanL_u,boundaryX,...
    preBCnum_nbrs, preBCthresh_mel, preBCthresh_xanL, preBCthresh_xanD);

cells_mel = cells_mel_u(logical_clean_mel,:);
cells_xanD = cells_xanD_u(logical_clean_xanD,:);
cells_xanL = cells_xanL_u(logical_clean_xanL,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Defining periodic boundary domain (in x) and computing pairwise
%%% cell-cell distances
periodicDistance = @(x,y)(sqrt((x(:,2)-y(:,2)).^2 + min((x(:,1)-y(:,1)).^2, (boundaryX -abs(x(:,1)-y(:,1))).^2)));
D_mel = pdist2(cells_mel, cells_mel, periodicDistance);        % distances between melanophores
D_xanD = pdist2(cells_xanD, cells_xanD, periodicDistance);     % distances between dense xanthophores
D_xanL = pdist2(cells_xanL, cells_xanL, periodicDistance);     % distances between loose xanthophores
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Saving distance data in .txt files for use when computing
%%% persistent homology.
if save_files
    writematrix(D_mel,sprintf('%s/%s_distances_mel_t=%d.txt',savefolder,simname,t),Delimiter=',');
    writematrix(D_xanD,sprintf('%s/%s_distances_xanD_t=%d.txt',savefolder,simname,t),Delimiter=',');
    writematrix(D_xanL,sprintf('%s/%s_distances_xanL_t=%d.txt',savefolder,simname,t),Delimiter=',');
end
pw_distances.M = D_mel;
pw_distances.Xd = D_xanD;
pw_distances.Xl = D_xanL;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end



