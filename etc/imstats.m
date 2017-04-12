function s = imstats(ref,x1)
%ref reference image
%x1 comparison image
%assumes images scaled to [0,1]
    s.MSE = (1/length(ref(:)))*norm(x1(:)-ref(:)).^2; 
    s.RMSE = sqrt(s.MSE);
    s.SNR = 20*log10(norm(ref(:))/norm(x1(:)-ref(:)));
    s.PSNR = 20*log10(norm(ref(:),Inf)/s.RMSE);  
    s.MSEabs = (1/length(ref(:)))*norm(abs(x1(:))-abs(ref(:))).^2; 
    s.RMSEabs = sqrt(s.MSEabs);      
    s.SNRabs = 20*log10(norm(ref(:))/norm(abs(x1(:))-abs(ref(:))));
    s.PSNRabs = 20*log10(norm(ref(:),Inf)/s.RMSEabs);    
    s.SSIM = ssim(im2uint8(abs(x1)),im2uint8(abs(ref)));
    s.HFEN = hfen(x1,ref);
    %s.relerr_inf = norm(x1(:)-x0(:),Inf)/norm(x0(:),Inf);
    %s.NMSE = (norm(x1(:)-x0(:)).^2)/(norm(x0(:)).^2);
    %s.NRMSE = sqrt(s.NMSE);
end