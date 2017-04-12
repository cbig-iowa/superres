name = 'phantom_brain';
% load data
load MR_brain_convert;

% set global parameters
params.hires = [256 256];
params.lores = [97 97];
params.filter_siz = [49 49];
params.cutoff = 1450;

% LSLP
params.LSLP.Niter = 500;
params.LSLP.lambda = Inf; 

% WTV
params.WTV.lambda = 1e-2;
params.WTV.Niter = 1000;

% TV
params.TV.lambda = 1e-5;
params.TV.Niter = 150;

% ncvxTV
params.ncvxTV.p = 0.5;
params.ncvxTV.lambda = 1e-4; 
params.ncvxTV.Niter = 800;
params.ncvxTV.beta = 100;
params.ncvxTV.betafac = 1.05;

% p = 0.5;
% lambda = 1e-3; %regularization parameter (typically in the range [1e-2,1], if original image scaled to [0,1])
% Niter = 200;  %number of iterations (typically 20-100 for Fourier inversion)
% beta = 100;
% betafac = 1.2;