%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Given a set of model parameters for the zebrafish model,
%compute persistence landscapes for all possible objectives.
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

function [output_file,psurfs,plands] = params_to_tda(sim_params,randomseed,T,tda_params,ts,savefolder,simname,saveextra)
%Simulate the Volkening zebrafish model, and run TDA analysis pipeline on the resulting pattern,
%which involves computing persistence homology using the Vietoris-Rips filtration,
%then compute persistence landscape, and optionaly persistence image.
% sim_params: struct containing parameter values for the zebrafish
%  simulation, pass in 0 or nan for default given by zebrafish_default_params()
% randomseed: rng seed for running the simulation
% T: endpoint of simulation (default in paper: T=45)
% tda_params: struct containing parameters for TDA analysis, pass in 0 for
%  default given by tda_default_params()
% ts: time points of interest to do TDA, can be [] (no TDA done) or scalar
%  (TDA done at one time), or array (TDA done at specified times).
%  Default is ts=1:T (do TDA at every time point)
% savefolder: folder to save results
% simname: name of simulation, used as prefix of all saved files
% saveextra: whether to save extra figures and intermediate files
%  (0: don't save, 1: save everything, 2: extra pattern pics only)
addpath('../zebrafish_abm');
fprintf('params_to_tda %s, begin: %s\n',simname,string(datetime('now'),'yyyy/MM/dd HH:mm:ss'));

% make sure savefolder have a slash at the end
savefolder = char(savefolder);
if savefolder(end) ~= '/'
    savefolder = strcat(savefolder,'/');
end

% default settings
if isstruct(sim_params)
    params=sim_params;
else
    params=zebrafish_default_params();
end
if ~isstruct(tda_params)
    tda_params=tda_default_params();
end

pland_kmax=tda_params.pland_kmax;
pland_numpts=tda_params.pland_numpts;
warning('off','MATLAB:Figure:FigureSavedToMATFile');
output_file = strcat(savefolder,simname,'.mat');

% simulation
tic;
zebrafish_model(params,T,savefolder,simname,randomseed);
fprintf('%s, Simualation took %.2f sec\n',simname,toc);

if saveextra>0
    for t=ts
        plot_zebrafish_pattern(savefolder,simname,t,4);
    end
else
    % always plot final pattern
    plot_zebrafish_pattern(savefolder,simname,T,4);
end

% TDA
num_ts=numel(ts);
formatted_data = cell(num_ts,1);
pw_distances = cell(num_ts,1);
barcodes=struct;
plands=struct;
psurfs=struct;
% How much domain in y-direction is cropped away for computing TDA
% default in Liu & Volkening, and Cleveland et al, is Ycrop=0.1
Ycrop=0.1;

for celltype=["Xd","M","Xl"]
    for dim =[0,1]
        for i=1:num_ts
            t=ts(i);
            fieldname=sprintf('%s_dim%d_t_%02d',celltype,dim,t);
            barcodes.(fieldname)=[];
            plands.(fieldname)=[];
            psurfs.(fieldname)=[];
        end
    end
end

for i=1:num_ts
    t=ts(i);

    % compute Vietoris-Rips persistence homology
    tic;
    [~,formatted_data{i},pw_distances{i}] = formatModelData(savefolder,simname,t,saveextra==1,Ycrop);
    fprintf('%s,t=%.2f, Formatting and distances took %.2f sec\n',simname,t,toc);
    pyenv(ExecutionMode="OutOfProcess");

    for celltype=["Xd","M","Xl"] % "M", "Xd", "Xl"
        for dim=[0,1]
            fieldname=sprintf('%s_dim%d_t_%02d',celltype,dim,t);
            if numel(pw_distances{i}.(celltype)) <= 1
                % there are 0 or 1 cells of this type, so no p-homology
                plands.(fieldname) = zeros(pland_kmax,pland_numpts);
                if do_pimg
                    psurfs.(fieldname) = zeros(tda_params.pimg_numpts, tda_params.pimg_numpts);
                end
                fprintf('%s, %s cells, t=%02d, dim=%d, empty p-homology\n',simname,celltype,t,dim);
                continue;
            end

            tic;
            [barcodes.(fieldname),pland_tmp_M]=pyrunfile("computePersistentHomology.py",{'barcodes','pland'},p_max=tda_params.pland_pmax,pland_step=tda_params.pland_numpts,dim=dim,D=py.numpy.array(pw_distances{i}.(celltype)));
            barcodes.(fieldname)=double(barcodes.(fieldname)); % convert from numpy to matlab type
            pland_tmp_M=double(pland_tmp_M);
            fprintf('%s, %s cells, t=%02d, dim=%d, Compute p-homology and landscape took %.2f sec\n',simname,celltype,t,dim,toc);

            % persistence landscape
            plands.(fieldname) = zeros(pland_kmax,pland_numpts);
            plands.(fieldname)(1:min(pland_kmax,size(pland_tmp_M,1)),1:pland_numpts) = pland_tmp_M(1:min(pland_kmax,size(pland_tmp_M,1)),1:pland_numpts);
            if saveextra==1
                fig_pland=plot_pland(plands.(fieldname),0);
                fig_barcode = plot_barcode(barcodes.(fieldname),0);
                saveas(fig_pland,sprintf('%s%s_pland_%s_t=%02d_dim%d.png',savefolder,simname,celltype,t,dim));
                saveas(fig_barcode,sprintf('%s%s_barcode_%s_t=%02d_dim%d.png',savefolder,simname,celltype,t,dim));
            end
        end
    end
end
clear('pw_distances'); % it takes too much memory
tic;
plands_file=sprintf('%s/%s_plands.mat',savefolder,simname);
save(plands_file,"-struct","plands");

barcode_file=sprintf('%s/%s_barcodes.mat',savefolder,simname);
save(barcode_file,"-struct","barcodes");
save(output_file,'-regexp', '^(?!(plands|psurfs|barcodes)$).','-append');% don't save plands and psurfs in the main file
fprintf('%s, saving took %.2f sec\n',simname,toc);
fprintf('params_to_tda %s, finished: %s\n',simname,string(datetime('now'),'yyyy/MM/dd HH:mm:ss'));
end
