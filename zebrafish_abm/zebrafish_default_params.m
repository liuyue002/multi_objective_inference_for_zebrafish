
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Default parameters for the zebrafish ABM.
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

function params = zebrafish_default_params()
%Default parameter values for zebrafish simulation

params.domain_reduce=1;
params.birthOption = 1; % 1 for wild-type, nacre, pfeffer, 3 for shady
params.cellsBornM = 350; % for wild-type, pfeffer, shady, 0 for nacre
params.cellsBornRandOrig = 500; % for wild-type, nacre, shady, 0 for pfeffer
params.cellsBornRandOrigM = 40; % for wild-type, shady, pfeffer, 0 for nacre

%%% Parameters related to cell movement - note the notation is Rxm = impact
%%% of xanthophores on melanophores, etc.
params.Rmm = 20;
params.Rxx = 16;
params.Rmx = 20;
params.Rxm = 15;
params.rmm = 50;
params.rxx = 45;
params.rxm = 67;
params.rmx = 67;
params.Axsnm = 0;
params.axsnm = 100;
params.Rxsnxsn = 16;
params.rxsnxsn = 60;
params.Rxxsn = 16;
params.rxxsn = 45;
params.Rxsnx = 16;
params.rxsnx = 45;
params.Aim = 0;
params.aim = 100;
params.Aixsn = 0;
params.Axi = 0;
params.axi = 20;
params.Ridid = 10;
params.Rilil = 32;
params.Ridil = 35;
params.Rilid = 32;
params.Aix = 5;
params.Rim = 5;
params.Rmi = 0;
params.ridid = 40;
params.rilil = 40;
params.ridil = 40;
params.rilid = 40;
params.aix = 90;
params.rim = 30;
params.rmi = 50;
params.Rilx = 0;
params.rilx = 60;
params.aixsn = params.aix;
params.Ailm = 0;
params.ailm = 20;
params.optionR = 50;

%%% Parameters related to cell birth as well as length scales used in cell
%%% interaction rules (birth, death, transitions in form)
params.localBernBM = 0.01;
params.MbornPropXtoMfarlb = 3.5; %beta
params.MbornAddPar = 3;    % alpha
params.MdeadPropMtoXfarlb = 2; % xi
params.MdeadPropXtoMloclb = 1.25; % mu
params.annulusOut = 250;
params.annulusIn = 210;
params.annulusOutI = 250;
params.annulusInI = 210;
params.ballR = 90;
params.dXM = 82;
params.ballR2 = 75;

%%% Overcrowding parameters to constrain random birth or bith at high
%%% densities
params.MoverPoplb = 4;
params.XoverPoplb = 6;
params.IoverPoplb = 7;
params.XsnoverPoplb = 1;
params.randomAmount = 0;

%%% Parameters used in iridophore and xanthophore transition rules
params.a = 2;
params.b = 1;
params.c = 3;
params.d = 9;
params.e = 3;
params.f = 3;
params.g = 5;
params.h = 2;

end