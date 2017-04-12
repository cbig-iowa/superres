function x_lo_nophase = phase_correction(x_lo,lores,hires,ind_lo)
%apply phase correction to low-resolution k-space samples
    if mod(lores(1),2)
        ramp1 = [(lores(1)+1)/2:-1:1, 1:(lores(1)-1)/2];
    else
        ramp1 = [lores(1)/2:-1:0, 1:(lores(1)/2-1)];
    end
    if mod(lores(2),2)
        ramp2 = [(lores(2)+1)/2:-1:1, 1:(lores(2)-1)/2];
    else
        ramp2 = [lores(2)/2:-1:0, 1:(lores(2)/2-1)];
    end
    x_pad = zeros(hires);
    x_pad(ind_lo) = x_lo;
    X = ifft2(x_pad);

    %x_filtered = ((ramp1.')*ramp2).*x_lo;
    x_filtered = x_lo;
    x_filtered_pad = zeros(hires);
    x_filtered_pad(ind_lo) = x_filtered;
    Xf = ifft2(x_filtered_pad);
    absXf = abs(Xf);
    absXf(absXf == 0) = 1;
    phase = exp(1j*angle(Xf));
    X_nophase = X.*conj(phase);
    x_nophase = fft2(X_nophase);
    x_lo_nophase = reshape(x_nophase(ind_lo),lores);
end