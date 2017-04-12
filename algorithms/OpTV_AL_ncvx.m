function [X, cost] = OpTV_AL_ncvx(b,lambda,A,At,res,Niter,p,beta,betafac)
% OPTV_AL: Solves TV regularized inverse problems with an
% ADMM scheme. Minimizes the cost function
% X* = argmin_X ||A(X)-b||_2^2 + lambda ||W |D(X)| ||_1
% where     X* = recovered image
%           A  = linear measurement operator
%           b  = (noisy) measurements
%           W  = diagonal weight matrix built from the edge mask
%           |D(X)| = gradient magnitude at each pixel
%
% Inputs:  A = function handle representing the forward
%               model/measurement operator
%          At = function handle representing the backwards model/
%               the transpose of the measurment operator.
%               (e.g. if A is a downsampling, At is a upsampling)                    
%          b =  a vector of measurements; should match the
%               dimensions of A(X)
%          lambda = regularization parameter that balances data fidelity
%               and smoothness. set lambda high for more smoothing.
%          siz = output image size, e.g. siz = [512,512]
%          Niter = is the number of iterations; should be ~500-1000
%         
% Output:  X = high-resolution output image
%          cost = array of cost function value vs. iteration

%beta = 10; %ADMM parameter (beta=10 should work well for images scaled between [0,1])

%Define AtA fourier mask
p_image = zeros(res); p_image(1,1) = 1;
AtA = fft2(At(A(p_image)));

%Define derivative operators D, Dt, and DtD
[D,Dt] = defDDt;
DtD = fft2(Dt(D(p_image)));

%Initialize X
Atb = At(b);
X = Atb;
DX = D(X);
G = zeros(size(DX));

% Begin alternating minimization alg.
cost = zeros(1,Niter);
for i=1:Niter
    %Shrinkage step
    Z = DX + G;
    Z1 = Z(:,:,1);
    Z2 = Z(:,:,2);
    AZ = sqrt(abs(Z1).^2 + abs(Z2).^2);
    shrinkZ = shrink_lp(AZ,1/beta,p); %shrinkage of gradient mag.
    Z(:,:,1) = shrinkZ.*Z1; 
    Z(:,:,2) = shrinkZ.*Z2;

    %Inversion step 
    F1 = fft2(2*Atb + lambda*beta*Dt(Z-G));
    F2 = lambda*beta*DtD + 2*AtA;
    X = ifft2(F1./F2);

    %Calculate error        
    DX = D(X);
    NDX = sqrt(abs(DX(:,:,1)).^2 + abs(DX(:,:,2)).^2);        
    diff = A(X)-b;
    cost(i) = norm(diff(:)).^2 + lambda*sum(NDX(:).^p); %objective function   
    %Lagrange multiplier update
    G = G + DX - Z;        
    beta = beta*betafac;
    beta = min(beta,1e20);
end
