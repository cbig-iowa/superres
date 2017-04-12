function y = proj(x,ind)
y = zeros(size(x));
y(ind) = x(ind);
y = y(:);
end