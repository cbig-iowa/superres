function edgemask = compute_sos_mask(V,s,filter_siz,outres,q,cutoff)
%% Build edgemask from annihilating filters U and weights s
kout = get_kspace_inds(outres);
ind_filter_out = get_lowpass_inds(kout,filter_siz);
mu2 = zeros(outres);
for j=cutoff:size(V,2)
    filter = zeros(outres);
    filter(ind_filter_out) = ifftshift(reshape(V(:,j),filter_siz));
    mu2 = mu2 + ((1/s(j))^q)*(abs(ifft2(filter)).^2);    
end
mu2 = mu2/max(abs(mu2(:)));
edgemask = sqrt(abs(mu2));
end