%% Load data
clc; clear all; close all;
addpath('./algorithms','./data','./etc','./exp');

%load data and parameters
%phantom_SL; %reproduce Figure 6--Shepp-Logan phantom
phantom_brain %reproduce Figure 6 -- brain phantom
%realdata %reproduce Figure 7 -- real data

X0 = ifft2(m); %input image
figure(1); imagesc(abs(X0)); colormap('gray'); colorbar; title('orig. mag.')
figure(2); imagesc(angle(X0)); colormap('gray'); colorbar; title('orig. phase')
%% Set target resolution
hires = params.hires;
lores = params.lores;
filter_siz = params.filter_siz;
%% Trim/pad original data
inres = size(m);
res = min(hires,inres);
k_hi = get_kspace_inds(hires);
k_in = get_kspace_inds(inres);
ind_hi = get_lowpass_inds(k_hi,res);
ind_orig = get_lowpass_inds(k_in,res);
x1 = zeros(hires);
x1(ind_hi) = m(ind_orig);
X1 = ifft2(x1); %hi-resolution
normfac = max(X1(:));
X1 = X1/normfac;
x1 = x1/normfac;
ind_lo = get_lowpass_inds(k_hi,lores);
x1_lo = reshape(x1(ind_lo),lores);
X1_lo = ifft2(x1_lo);
x_lo = x1_lo;
x_ref_lo = x1_lo;
x_ref = x1;
X_ref = ifft2(x1);
figure(1); imagesc(abs(X1_lo)); colormap('gray'); colorbar; title('low-res mag.')
figure(2); imagesc(angle(X1_lo)); colormap('gray'); colorbar; title('low-res phase')
%% Add noise
if isfield(params,'add_noise')
    sig = params.add_noise.sig;
    x1_lo = x_ref_lo + sig*(randn(size(x_lo)) + 1i*randn(size(x_lo)));
    noiseSNR = 20*log10(norm(x_ref_lo)/norm(x_ref_lo-x1_lo));
    fprintf('Sample SNR = %5.2f dB\n',noiseSNR);
    
    figure(1); imagesc(abs(ifft2(x1_lo))); colormap('gray'); colorbar; title('noisy low-res mag.')
    figure(2); imagesc(angle(ifft2(x1_lo))); colormap('gray'); colorbar; title('noisy low-res phase')    
end
%% Optional phase correction for real MR data
if isfield(params,'phase_correct')
    if params.phase_correct
        x_lo = phase_correction(x1_lo,lores,hires,ind_lo);
        k_hi2 = get_kspace_inds(2*hires);
        ind_hi2 = get_lowpass_inds(k_hi2,hires);
        x_ref = phase_correction(x1,hires,2*hires,ind_hi2);
        %X_ref = ifft2(x_ref);
        X_ref = abs(X1);
    end
    X_lo = ifft2(x_lo);  
    figure(3); imagesc(abs(X_lo)); colormap('gray'); colorbar; title('phase corrected low-res mag.')
    figure(4); imagesc(angle(X_lo)); colormap('gray'); colorbar; title('phase corrected low-res phase')
    figure(13); imagesc(abs(X_ref)); colormap('gray'); colorbar; title('phase corrected high-res mag.')
    figure(14); imagesc(angle(X_ref)); colormap('gray'); colorbar; title('phase corrected high-res phase')          
end
%norm(x_lo(:)-x1_lo(:))/norm(x1_lo(:))
%figure; imagesc((abs(X_lo)-abs(X1_lo))); colorbar;
%% View zero-padded IFFT reconstruction for reference
sampmask = zeros(res);
sampmask(ind_lo) = 1;
x_ifft = zeros(hires);
x_ifft(ind_lo) = x1_lo;
X_IFFT = ifft2(x_ifft);
figure(5); imagesc(abs(X_ref),[0 1]); colorbar; colormap('gray'); title('orig. at target resolution')
figure(6); imagesc(abs(X_IFFT),[0 1]); colorbar; colormap('gray'); title('zero-padded IFFT low resolution')
figure(7); imagesc(abs(fftshift(x_ifft)).^(1/4)); colorbar; colormap('jet'); title('zero-padded k-space data');
figure(8); imagesc(abs(abs(X_IFFT)-abs(X_ref)),[0,0.2]); colorbar; colormap('jet'); title('IFFT, error image (magnitude only)');
stats.IFFT = imstats(X_ref,X_IFFT);
fprintf('IFFT results:\t SNR=%2.2f dB, HFEN=%1.4f, SSIM=%1.4f\n',stats.IFFT.SNRabs,stats.IFFT.HFEN,stats.IFFT.SSIM);
%% STAGE 1: Compute annihilating subspace/sos polynomial

%% Build Toeplitz lifting operator T
k_lo = get_kspace_inds(lores);
dz = zeros([lores,2]);
dz(:,:,1) = reshape(1j*2*pi*k_lo(1,:),lores);
dz(:,:,2) = reshape(1j*2*pi*k_lo(2,:),lores);
[TL, ThL] = define_toeplitz_operators(dz,lores,filter_siz,0);
%% Compute SVD of lifted data
[U,S,V] = svd(TL(x_lo),'econ');
s = diag(S);
figure(8); plot(s/s(1)); title('normalized singular values of T(fhat)');
%% set rank cutoff to identify annihilating subspace
params.cutoff = floor(0.5*prod(params.filter_siz));
cutoff = params.cutoff;
sos = compute_sos_mask(V,s,filter_siz,hires,0,cutoff);
figure(2); imagesc(sos,[0,1]); colorbar; colormap('default'); title('sum of squares mask');
VN = V(:,cutoff:end); %annhilating subspace
%% STAGE 1.5 (optional): Cadzow denoise to refine mask
% default parameters
% params.cadzow.lambda = 5e-1;
% params.cadzow.iter = 5;
% params.cutoff = floor(0.5*prod(params.filter_siz));
if isfield(params,'cadzow')
lambda = params.cadzow.lambda;
r = params.cutoff; %rank threshold

x = x_lo;
Tx = TL(x);
Z = zeros(size(Tx)); %aux variable
ThT = ThL(TL(ones(lores)));

cost = [];
constraint = [];
iter = params.cadzow.iter;
for i=1:iter
    Tx = TL(x);
    [U,S,V] = svd(Tx,'econ');
    s = diag(S);
    
    figure(4); plot(log10(s)); title('singular values of T(x), log scale');
    sos = compute_sos_mask(V,s,filter_siz,[256,256],0,r);
    figure(2); imagesc(sos,[0,1]); colorbar; colormap('default'); title('sum of squares mask'); drawnow;
    
    ranknorm = 0.5*norm(s((r+1):end))^2;    
    s((r+1):end) = 0;
    
    Z = U*diag(s)*V';    
    x = (x_lo + lambda*ThL(Z))./(1+lambda*ThT);
    
    diff = x-x_lo;
    thiscost = 0.5*norm(diff(:)).^2 + lambda*ranknorm;
    cost = [cost,thiscost];        

    figure(6); imagesc(abs(ifft2(x))); colorbar; drawnow;
    figure(7); plot(cost); title('cost'); drawnow;
end

Tx = TL(x);
[U,S,V] = svd(Tx,'econ');
s = diag(S);
x_lo = x;

cutoff = params.cutoff;
sos = compute_sos_mask(V,s,filter_siz,hires,0,cutoff);
figure(2); imagesc(sos,[0,1]); colorbar; colormap('default'); title('sum of squares mask'); drawnow
VN = V(:,cutoff:end); %annhilating subspace
end
%% STAGE 2: Recover super-resolved image

%% Method 1: Fourier domain recovery via least squares linear prediction (LSLP)
%Solves linear least squares problem 
% min_x ||sos \ast x||_2^2 + lambda*||P(x)-b||_2^2 with LSQR

% default parameters
%params.LSLP.Niter = 500;
%params.LSLP.lambda = 1e-1;
Niter = params.LSLP.Niter;   %number of LSQR iterations, typically 200-500 is sufficient
lambda = params.LSLP.lambda; %Set lambda=Inf to solve with equality constraint on data.
                             %Otherwise, for noisy data, typically in range [1e-1,100];
                             
b = x_lo; %cadzow denoised recovery
[x_lslp,flag,relres,iter,resvec] = OpLSLP(b,hires,lores,VN,filter_siz,Niter,lambda);
X_LSLP = ifft2(x_lslp);
figure(100); imagesc(abs(X_LSLP),[0,1]); colorbar; colormap('gray'); title('LSLP solution');
figure(110); imagesc(fftshift(log(1+abs(x_lslp)))); colormap('jet'); colorbar; title('LSLP cartoon solution, k-space'); 
figure(120); imagesc(abs(abs(X_LSLP)-abs(X_ref)),[0,0.2]); colormap('default'); colorbar; title('LSLP, error image (magnitude only)');
stats.LSLP = imstats(X_ref,X_LSLP);
fprintf('LSLP results:\t SNR=%2.2f dB, HFEN=%1.4f, SSIM=%1.4f, (lambda=%2.3f)\n',stats.LSLP.SNRabs,stats.LSLP.HFEN,stats.LSLP.SSIM,lambda);
%% Method 2: Spatial domain recovery via weighted total variation (WTV)
% Solved the convex weighted total variation minimization problem:
% min_x ||W.*Dx||_1 + lambda*||P(x)-b||_2^2 with a primal-dual algorithm
% where W = weights determined by sos polynomial,
% Dx = discrete gradient

% default parameters
%params.WTV.lambda = 1e-3;
%params.WTV.Niter = 1000;
[A,At] = defAAt_fourier(ind_lo, hires); %Define function handles for fourier projection operators                                                
b = (1/sqrt(prod(hires)))*x1_lo(:).'; %low-resolution fourier samples
lambda = params.WTV.lambda; %regularization parameter (typically in the range [1e-2,1], if original image scaled to [0,1])
Niter = params.WTV.Niter;  %number of iterations (typically 500-1000 for Fourier inversion)
[X_WTV, cost] = OpWeightedTV_PD(b,sos,lambda,A,At,hires,Niter); %see comments inside OpWeightedTV_PD.m
figure(20); imagesc(abs(X_WTV),[0,1]); colorbar; colormap('gray'); title('WTV solution');
figure(21); imagesc(fftshift(log(1+abs(fft2(X_WTV))))); colorbar; title('WTV solution, k-space'); 
figure(22); imagesc(abs(abs(X_WTV)-abs(X_ref)),[0,0.2]); colorbar; title('WTV, error image (magnitude only)');
stats.WTV = imstats(X_ref,X_WTV);
fprintf('WTV results:\t SNR=%2.2f dB, HFEN=%1.4f, SSIM=%1.4f\n',stats.WTV.SNRabs,stats.WTV.HFEN,stats.WTV.SSIM);
%% Comparisons with other methods

%% Run Standard TV algorithm (ADMM implementation)
% default parameters
%params.TV.lambda = 1e-2;
%params.TV.Niter = 200;
[A,At] = defAAt_fourier(ind_lo, hires); %Define function handles for fourier projection operators                                                
b = (1/sqrt(prod(hires)))*x1_lo(:).';%A(X1); %low-resolution fourier samples
lambda = params.TV.lambda; %regularization parameter (optimal range typically [1e-2,1], if original image scaled to [0,1])
Niter = params.TV.Niter; %number of iterations (typically 20-100 for Fourier inversion)
[X_TV, cost] = OpTV_AL(b,lambda,A,At,hires,Niter); %see comments inside OpWeightedTV_PD.m
figure(30); imagesc(abs(X_TV),[0,1]); colormap('gray'); colorbar; title('TV solution');
figure(31); imagesc(fftshift(log(1+abs(fft2(X_TV))))); colorbar; title('TV solution, k-space'); 
figure(32); imagesc(abs(abs(X_TV)-abs(X_ref)),[0,0.2]); colorbar; title('TV, error image (magnitude only)');
stats.TV = imstats(X_ref,X_TV);
fprintf('TV results:\t SNR=%2.2f dB, HFEN=%1.4f, SSIM=%1.4f\n',stats.TV.SNRabs,stats.TV.HFEN,stats.TV.SSIM);
%% Run Non-convex TV algorithm (Chartrand l_p penalty -- ADMM implementation)
% default parameters
%params.ncvxTV.p = 0.5;
%params.ncvxTV.lambda = 1e-3; 
%params.ncvxTV.Niter = 120;
%params.ncvxTV.beta = 100;
%params.ncvxTV.betafac = 1.3;
[A,At] = defAAt_fourier(ind_lo, hires); %Define function handles for fourier projection operators                                                
b = (1/sqrt(prod(hires)))*x1_lo(:).'; %low-resolution fourier samples
p = params.ncvxTV.p; %use non-convex l^p norm p < 1. Typically p in range [0.25,0.5] works well.
Niter = params.ncvxTV.Niter;  %number of iterations (typically 20-100 for Fourier inversion)
beta = params.ncvxTV.beta;  %ADMM parameter
betafac = params.ncvxTV.betafac; %continuation factor on ADMM parameter (in range [1.1,2]);
lambda = params.ncvxTV.lambda; %regularization parameter (optimal range typically [1e-5,1e-1], if original image scaled to [0,1])
[X_ncvxTV, cost] = OpTV_AL_ncvx(b,lambda,A,At,hires,Niter,p,beta,betafac); %see comments inside OpWeightedTV_PD.m
figure(40); imagesc(abs(X_ncvxTV),[0,1]); colormap('gray'); colorbar; title('nonconvex TV soution');
figure(41); imagesc(fftshift(log(1+abs(fft2(X_ncvxTV))))); colorbar; title('nonconvex TV solution, k-space'); 
figure(42); imagesc(abs(abs(X_ncvxTV)-abs(X_ref)),[0,0.2]); colorbar; title('nonconvex TV, error image (magnitude only)');
figure(43); plot(cost); xlabel('iteration'); ylabel('cost'); %non-convex alg--ensure cost is decreasing
stats.ncvxTV = imstats(X_ref,X_ncvxTV);
fprintf('ncvxTV results:\t SNR=%2.2f dB, HFEN=%1.4f, SSIM=%1.4f\n',stats.ncvxTV.SNRabs,stats.ncvxTV.HFEN,stats.ncvxTV.SSIM);