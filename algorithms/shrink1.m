function Y = shrink1(X,tau)
%Returns shrunk X using l^1 shrinkage
Z = abs(X);
M = ones(size(X));
M(Z==0) = 0;
Z(Z==0) = 1;
S = max(1 - tau./Z, 0).*M;
Y = S.*X;

end


