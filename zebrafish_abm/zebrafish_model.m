function [output_file] = zebrafish_model(params,T,savefolder,simname,randomseed)
%Simulate the zebrafish model.
%The code is adapted from Volkening & Sandstede (2018), and provided
%separately upon request. 
% params: values for model parameters. Defaults given in
% zebrafish_default_params.m
% T: end point of simulation. Default: 45 (=65 days post fertilisation)
% simname: a prefix in the name of all output files
% randomseed: rng seed for the simulation
output_file = '';
end