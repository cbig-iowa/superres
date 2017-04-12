function A = lsqrA(sos,dz,overres,ind_samples_comp)
%FULLR Summary of this function goes here
    ndz = size(dz,3);
    sosrep = repmat(sos,[1,1,ndz]);
    
    function y = afun(x,transp_flag)
       if strcmp(transp_flag,'transp')      % y = A'*x
          z = sum(conj(dz).*fft2(sosrep.*reshape(x,[overres,ndz])),3);
          y = z(ind_samples_comp);
          y = y(:)/sqrt(prod(overres));
       elseif strcmp(transp_flag,'notransp') % y = A*x
          z = zeros(overres);
          z(ind_samples_comp) = x;          
          y = sosrep.*ifft2(repmat(z,[1,1,ndz]).*dz);
          y = sqrt(prod(overres))*y(:);
       end
    end

    A = @(x,transp_flag) afun(x,transp_flag);
end

