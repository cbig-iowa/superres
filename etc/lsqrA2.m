function A = lsqrA2(sos,dz,overres,ind_samples,lambda)
%FULLR Summary of this function goes here
    ndz = size(dz,3);
    sosrep = repmat(sos,[1,1,ndz]);
    N = prod([overres,ndz]);
    
    function y = afun(x,transp_flag)
       if strcmp(transp_flag,'transp')      % y = A'*x
          x1 = x(1:N);
          x2 = x((N+1):end);          
          y1 = sum(conj(dz).*fft2(sosrep.*reshape(x1,[overres,ndz])),3)/sqrt(prod(overres));
          y2 = zeros(overres);
          y2(ind_samples) = sqrt(lambda)*x2;
          y = y1+y2;
          y = y(:);
       elseif strcmp(transp_flag,'notransp') % y = A*x         
          y1 = sosrep.*ifft2(repmat(reshape(x,overres),[1,1,ndz]).*dz);
          y1 = sqrt(prod(overres))*y1(:);
          y2 = sqrt(lambda)*x(ind_samples);
          y = [y1(:); y2(:)];
       end
    end

    A = @(x,transp_flag) afun(x,transp_flag);
end

