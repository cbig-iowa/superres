function [X, cost] = OpWeightedTV_PD(b,edgemask,lambda,A,At,res,Niter)
% OPWEIGHTEDTV_PD: Solves weighted TV regularized inverse problems with a
% primal-dual scheme. Minimizes the cost function
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
%  
% Note this implementation is an adaptation of the algorithm in the paper:
%   Zhu, M., & Chan, T. (2008). "An efficient primal-dual hybrid 
%       gradient algorithm for total variation image restoration." 
%       UCLA CAM Report, 08-34.

%Define AtA fourier mask
p_image = zeros(res,'double'); p_image(1,1) = 1;
AtA = fft2(At(A(p_image)));

%Define derivative operator
[D,Dt] = defDDt;
%Defined weighted derivative operators
Wbig = repmat(edgemask,[1,1,2]);
WD = @(x) Wbig.*D(x);
WDt = @(x) Dt(Wbig.*x);

%Initialize variables
Atb = At(b);
X = Atb;
WDX = WD(X);
P = zeros(size(WDX));

tau = 0.2;   %gradient step-sizes--leave these fixed
theta = 0.5;

lambda2 = 2/lambda;

% Begin alternating minimization alg.
cost = zeros(1,Niter);
for i=1:Niter    
    %Dual Step
    P = projInfty(P + tau*WDX);        

    %Primal Step
    X = ifft2(fft2(X - theta*WDt(P) + theta*lambda2*Atb)./(1 + theta*lambda2*AtA));
    
    %Calculate cost function      
    WDX = WD(X);
    NWDX = sqrt(abs(WDX(:,:,1)).^2 + abs(WDX(:,:,2)).^2);
    diff = A(X)-b;
    cost(i) = norm(diff(:)).^2 + lambda*sum(NWDX(:)); %cost function

    %Update gradient steps
    tau = tau + 0.08;
    theta = 0.5/tau;
end
