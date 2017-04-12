function y = projunitball(x,~)

normx = sqrt(abs(x(:,:,1)).^2 + abs(x(:,:,2)).^2);
normx(normx <= 1) = 1;
normx = repmat(normx,[1,1,2]);
y = x./normx;

end