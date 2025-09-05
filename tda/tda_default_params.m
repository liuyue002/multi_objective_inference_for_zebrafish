function [tda_params] = tda_default_params()
%Default parameters for TDA
tda_params.pland_pmax=600;   % maximum persistence to consider for computing p-landscape
tda_params.pland_numpts=600; % number of grid points for computing p-landscape
tda_params.pland_kmax=300;   % how many degree of p-landscape fto compute
tda_params.pimg_numpts=301;  % number of grid points for computing p-image
tda_params.pimg_bmax=600;    % max birth time for p-image
tda_params.pimg_pmax=600;    % max persistence for p-image

% smoothing function and weights for p-image
sig = 20; % smoothing
tda_params.phi = @(x,y,u,v) (1/(2*pi*sig^2))*exp(-((x-u)^2+(y-v)^2)/(2*sig^2)); % kernel
tda_params.w = @(b,p) 1/(1+exp(-(p-200)/25));
% precompute the convolution mask from phi (this saves time if doing p-img many times)
xx2=linspace(-tda_params.pimg_bmax,tda_params.pimg_bmax,tda_params.pimg_numpts*2-1);
yy2=linspace(-tda_params.pimg_pmax,tda_params.pimg_pmax,tda_params.pimg_numpts*2-1);
[X2,Y2] = meshgrid(xx2,yy2);
% assumes phi is symmeric
phi0=@(x,y) tda_params.phi(x,y,0,0);
tda_params.phimask = arrayfun(phi0,X2,Y2);
end