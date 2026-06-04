%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Wrapper function for running simulation and TDA computation.
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

function iridophore_sim(irid_params,mutant,randomseed,savefolder)

sim_params = zebrafish_default_params();
sim_params.c = irid_params(1);
sim_params.d = irid_params(2);
sim_params.e = irid_params(3);
sim_params.f = irid_params(4);
sim_params.g = irid_params(5);
sim_params.h = irid_params(6);

switch mutant
    case "wt"
        % do nothing
    case "pfeffer"
        sim_params.cellsBornRandOrig = 0;
    case "nacre"
        sim_params.cellsBornM = 0;
        sim_params.cellsBornRandOrigM = 0;
    case "shady"
        sim_params.birthOption = 3;
    otherwise
        error("No such mutant type");
end

T=45;
tda_params=tda_default_params();
ts=1:45;
simname=sprintf('zebrafish_iridophore_[%02d,%02d,%02d,%02d,%02d,%02d]_%s,rng=%d',irid_params,mutant,randomseed);
[~,~,~] = params_to_tda(sim_params,randomseed,T,tda_params,ts,savefolder,simname,0,1);
end
