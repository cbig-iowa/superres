name = 'realdata_cor2_cc_noisy4';
% load data
load MR_realbrain_cor2_cc;

% set global parameters
params.hires = [256,256];
params.lores = [128,128];
params.filter_siz = [64,64];
params.cutoff = floor(0.5*prod(params.filter_siz));
params.phase_correct = true;

params.add_noise.sig = 10;

% Cadzow denosing
params.cadzow.lambda = 5e-1;
params.cadzow.iter = 5;

% LSLP
params.LSLP.Niter = 500;
params.LSLP.lambda = 2; 

% WTV
params.WTV.lambda = 1e-5;
params.WTV.Niter = 1000;

% TV
params.TV.lambda = 1e-5;
params.TV.Niter = 200;

% ncvxTV
params.ncvxTV.p = 0.5;
params.ncvxTV.lambda = 1e-5; 
params.ncvxTV.Niter = 150;
params.ncvxTV.beta = 100;
params.ncvxTV.betafac = 1.3;