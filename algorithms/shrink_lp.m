function Y = shrink_lp(X,tau,p)
%Returns thresholded magnitudes of X.
Z = abs(X);
M = ones(size(X));
M(Z==0) = 0;
Z(Z==0) = 1;
Y = max(1 - tau./(Z.^(2-p)), 0).*M; %lp, p < 1 shrinkage formula

end


