name = 'phantom_SL';
% load data
load MR_SL_convert;

% set global parameters
params.hires = [256 256];
params.lores = [49 65];
params.filter_siz = [25 33];
params.cutoff = 330;

% LSLP
params.LSLP.Niter = 300;
params.LSLP.lambda = Inf; 

% WTV
params.WTV.lambda = 1e-3;
params.WTV.Niter = 1000;

% TV
params.TV.lambda = 1e-5;
params.TV.Niter = 150;

% ncvxTV
params.ncvxTV.p = 0.25;
params.ncvxTV.lambda = 1e-5; 
params.ncvxTV.Niter = 800;
params.ncvxTV.beta = 100;
params.ncvxTV.betafac = 1.05;