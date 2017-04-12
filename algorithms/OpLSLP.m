function [x,flag,relres,iter,resvec] = OpLSLP(b,hires,lores,V,filter_siz,Niter,lambda)
% OPLSLP Solves Least-squares linear prediction problem
% min_{x} ||sos \ast x||^2_2 + lambda*||P(x) - b||_2^2
% where x is extrapolated fourier data
% sos is the sum-of-squares annihilating filter built from annihliating
% subspace 
% P is projection onto the sampling set
% b are the known low-resolution fourier samples
% and lambda is a regularization parameter
% setting lambda = Inf solves equality constrained version:
% min_{ghat} ||sos \ast ghat||^2_2  s.t. P(ghat)=b
%
% Inputs: b = low-resolution data
%         hires = hi resolution dimensions
%         lores = lo resolution dimnesions
%         V = annihilating subspace -- computed from SVD of T(fhat)
%         filter_siz = size of annhilating filter
%         Niter = # of lsqr iterations (typically 200-400 works well)
%         lambda = regularization parameter (note setting lambda=Inf solves
%         equaltity constrained problem)
% Outputs: x = extrapolated fourier data on hi-resolution grid
%          flag,...,resvec = outputs from lsqr algorithm

overres = 2*hires; %set "over-resolved" grid at which to solve the problem--2*hires is a good rule of thumb
k_over = get_kspace_inds(overres);
ind_full = get_lowpass_inds(k_over,hires);
ind_lp   = get_lowpass_inds(k_over,lores);
mask_pad = zeros(overres);
mask_pad(ind_lp) = 1;
x0_pad = zeros(overres);
x0_pad(ind_lp) = b;
ind_samples_comp = find(mask_pad==0);
dz = zeros([overres,2]);
dz(:,:,1) = reshape(1j*2*pi*k_over(1,:),overres)/sqrt(prod(overres));
dz(:,:,2) = reshape(1j*2*pi*k_over(2,:),overres)/sqrt(prod(overres));
sos_pad = compute_sos_mask(V,ones(1,size(V,2)),filter_siz,overres,0,1);

if lambda == Inf %solve equality constraintd
    A = lsqrA(sos_pad,dz,overres,ind_samples_comp);
    B = -repmat(sos_pad,[1,1,size(dz,3)]).*(ifft2(repmat(x0_pad,[1,1,size(dz,3)]).*dz));
    B = sqrt(prod(overres))*B(:);

    pre_mask = sqrt(sum(k_over(:,ind_samples_comp).^2)).'; %precondtioner
    PRE = @(x,transp_flag) x./pre_mask;
    %
    [y,flag,relres,iter,resvec] = lsqr(A,B,1e-12,Niter,PRE);
    y_pad = zeros(overres);
    y_pad(ind_samples_comp) = y;
    x_pad = y_pad + x0_pad;
else
    A2 = lsqrA2(sos_pad,dz,overres,ind_lp,lambda);
    B1 = zeros(size(dz));
    B2 = sqrt(lambda)*b;
    B = [B1(:);B2(:)];

    pre_mask = sqrt(sum(k_over.^2)).';
    pre_mask(1) = 1;
    PRE2 = @(x,transp_flag) x./pre_mask; %precondtioner
    %
    [x_pad,flag,relres,iter,resvec] = lsqr(A2,B,1e-12,Niter,PRE2);
end

x = reshape(x_pad(ind_full),hires); %truncate over-resolved data to hi-resolution grid
end