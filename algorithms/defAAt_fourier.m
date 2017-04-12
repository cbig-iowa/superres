%Defines Fourier Undersampling operators A, At for Compressed Sensing MRI
%Inputs: ind_samples = indicies of sampling locations, res = resolution
%Outputs: A = measurement operator, At = measurement operator transpose
function [A,At] = defAAt_fourier(ind_samples,res)    
    A = @(z) A_fhp(z,ind_samples,res);
    At = @(z) At_fhp(z,ind_samples,res);

    function out = A_fhp(z, ind_samples, res)
        p = 1/sqrt(res(1)*res(2))*fft2(z);
        out = p(ind_samples);
    end

    function out = At_fhp(z,ind_samples, res)
        p = zeros(res,'double');
        p(ind_samples) = z;
        out = sqrt(res(1)*res(2))*ifft2(p);
    end
end