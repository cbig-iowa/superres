%% Code to generate figure 3
%% pwc image 1: non-intersecting regions
%Uncomment code to re-generate minimal polynomial and PWC function

% K0 = 1; %define a 3x3 filter
% m = 4*2048; %spatial grid size
% chi = zeros(m);
% xval = [ -1 0 1; -1 0 1; -1 0 1];
% yval = xval.';
% trans = @(x,y) exp(-1i*pi*xval*x).*exp(1i*pi*yval*y);
% 
% muhat = zeros(2*K0+1,2*K0+1,2*K0+1);
% muhat(:,:,1) = trans(0.3,0).*[0 1i 0; -1i -3 1i; 0 -1i 0];   %blob1
% muhat(:,:,2) = trans(1,0.5).*[-1 2i 1; -2i -5.0 2i; -1 -2i -1]; %blob2
% muhat(:,:,3) = trans(0,0.9).*[0.5 1i 0.5; -2i -4.5 2i; 0.5 -1i 0.5]; %blob3
% 
% amps = [0.2, -0.3, 0.8]; %region amplitudes
% realc = 1;
% for i=1:3
% muhat0 = muhat(:,:,i);
% muhat0 = fftshift(fft2(real(ifft2(ifftshift(muhat0))))); %make polynomial real-valued
% 
% ind = mod((1:(2*K0+1))-(K0+1),m)+1;
% J = zeros(m);
% J(ind,ind) = muhat0;
% I = m^2*real(ifft2(J));
% chi = chi + amps(i)*(I>0); 
% realc = conv2(realc,muhat0);
% end
% chi = flipud(fftshift(chi));
% 
% %compute fourier data
% K = 10; %square sampling window from -2K:2K by -2K:2K
% fhat0 = fft2(chi);
% ind = mod((1:(4*K+1))-(2*K+1),m)+1;
% fhat = fhat0(ind,ind);
% kx = -(2*K):(2*K);
% ky = kx.';
% 
% fxhat = fhat.*repmat(1i*2*pi*kx,[4*K+1,1]); %low-pass fourier coefficients of \partial_x f 
% fyhat = fhat.*repmat(1i*2*pi*ky,[1,4*K+1]); %low-pass fourier coefficients of \partial_x f 
% 
% fsmall = imresize(chi,[512,512],'nearest');
% save('pwc1.mat','realc','fxhat','fyhat','fsmall');
% clear I J chi;
%% pwc image 2: intersecting regions
% K0 = 1;
% m = 4*2048;
% chi = zeros(m);
% %muhat = randn(2*K+1,2*K+1) + 1i*randn(2*K+1,2*K+1);
% xval = [ -1 0 1; -1 0 1; -1 0 1];
% yval = xval.';
% trans = @(x,y) exp(-1i*pi*xval*x).*exp(1i*pi*yval*y);
% 
% muhat = zeros(3,3,3);
% muhat(:,:,1) = trans(0.5,0.1).*[0 1i 0; -1i -1.5 1i; 0 -1i 0];
% muhat(:,:,2) = trans(0.75,0.5).*[-1 2i 1; -2i -3.0 2i; -1 -2i -1];
% muhat(:,:,3) = trans(0.2,0.6).*[0.5 1i 0.5; -2i -2.5 2i; 0.5 -1i 0.5];
% 
% amps = [0.2, -0.3, 0.8];
% realc = 1;
% for i=1:3
% muhat0 = muhat(:,:,i);
% muhat0 = fftshift(fft2(real(ifft2(ifftshift(muhat0)))));
% 
% ind = mod((1:(2*K0+1))-(K0+1),m)+1;
% J = zeros(m);
% J(ind,ind) = muhat0;
% I = m^2*real(ifft2(J));
% chi = chi + amps(i)*(I>0);
% realc = conv2(realc,muhat0);
% end
% chi = flipud(fftshift(chi));
% a = unique(chi);
% a = a(a~=0);
% a = sort(a);
% b = [0.2 -0.1 -0.5 0.6 0.9 -0.2 0.4];
% for i=1:length(a)
%         chi(chi==a(i)) = b(i);
% end
% figure(1); imshow(chi,[]);
% chi = round(chi*10)/10;
% 
% %compute fourier data
% K = 10; %square sampling window from -2K:2K by -2K:2K
% fhat0 = fft2(chi);
% ind = mod((1:(4*K+1))-(2*K+1),m)+1;
% fhat = fhat0(ind,ind);
% kx = -(2*K):(2*K);
% ky = kx.';
% 
% fxhat = fhat.*repmat(1i*2*pi*kx,[4*K+1,1]); %low-pass fourier coefficients of \partial_x f 
% fyhat = fhat.*repmat(1i*2*pi*ky,[1,4*K+1]); %low-pass fourier coefficients of \partial_x f 
% 
% fsmall = imresize(chi,[512,512],'nearest');
% save('pwc2.mat','realc','fxhat','fyhat','fsmall');
% clear I J chi;
%% pwc image 3: intersecting regions
% K0 = 1;
% m = 4*2048;
% chi = zeros(m);
% %muhat = randn(2*K+1,2*K+1) + 1i*randn(2*K+1,2*K+1);
% xval = [ -1 0 1; -1 0 1; -1 0 1];
% yval = xval.';
% trans = @(x,y) exp(-1i*pi*xval*x).*exp(1i*pi*yval*y);
% 
% muhat = zeros(3,3,2);
% %muhat(:,:,1) = trans(0.5,0.1).*[0 1i 0; -1i -1.5 1i; 0 -1i 0];
% %muhat(:,:,1) = trans(0.75,0.5).*[-1 2i 1; -2i -3.0 2i; -1 -2i -1];
% muhat(:,:,1) = trans(0.5,0.5).*[0.1 1i 0.1; -2i -2.5 2i; 0.1 -1i 0.1];
% muhat(:,:,2) = trans(-0.5,0.5).*rot90([0.1 1i 0.1; -2i -2.5 2i; 0.1 -1i 0.1],1);
% 
% amps = [1, 2];
% realc = 1;
% for i=1:2
% muhat0 = muhat(:,:,i);
% muhat0 = fftshift(fft2(real(ifft2(ifftshift(muhat0)))));
% 
% ind = mod((1:(2*K0+1))-(K0+1),m)+1;
% J = zeros(m);
% J(ind,ind) = muhat0;
% I = m^2*real(ifft2(J));
% chi = chi + amps(i)*(I>0);
% realc = conv2(realc,muhat0);
% end
% 
% chi = flipud(fftshift(chi));
% 
% b = [-0.4147 - 0.6249i, 1.9122 - 1.1687i,  -0.3909 + 0.3926i, 0.4092 + 1.3018i, -1.1424 - 0.5936i];
% chi((chi==3)) = b(1);
% mask1 = chi==1;
% mask = mask1;
% mask(1:floor(m/2),:) = 0;
% chi(mask) = b(2);
% mask = mask1;
% mask(floor(m/2):end,:) = 0;
% chi(mask) = b(3);
% mask1 = chi==2;
% mask = mask1;
% mask(:,1:floor(m/2)) = 0;
% chi(mask) = b(4);
% mask = mask1;
% mask(:,floor(m/2):end) = 0;
% chi(mask) = b(5);
% 
% figure(1); imshow(abs(chi),[]);
% chi = round(chi*10)/10;
% 
% %compute fourier data
% K = 10; %square sampling window from -2K:2K by -2K:2K
% fhat0 = fft2(chi);
% ind = mod((1:(4*K+1))-(2*K+1),m)+1;
% fhat = fhat0(ind,ind);
% kx = -(2*K):(2*K);
% ky = kx.';
% 
% fxhat = fhat.*repmat(1i*2*pi*kx,[4*K+1,1]); %low-pass fourier coefficients of \partial_x f 
% fyhat = fhat.*repmat(1i*2*pi*ky,[1,4*K+1]); %low-pass fourier coefficients of \partial_x f 
% 
% fsmall = imresize(abs(chi),[512,512],'nearest');
% save('pwc3.mat','realc','fxhat','fyhat','fsmall');
% clear I J chi;
%% Load pre-computed data
load pwc1;
%load pwc2;
%% Plot piecewise constant function f and the minimal polynomial mu
K0 = 3;
ind = fliplr(mod((1:(2*K0+1))-(K0+1),512)+1);
muplot = zeros(512);
muplot(ind,ind) = realc;
muplot = fliplr(fftshift(real(ifft2(muplot))));
muplot = muplot/max(abs(muplot(:)));
figure(1); imagesc(fsmall,[-0.5,1]); colormap parula; colorbar; axis image; axis off; title('original PWC image');
figure(2); imagesc(muplot,[-0.15,0.15]); colormap parula; colorbar; axis image; axis off; title('original minimal polynomial');
hold on; contour(muplot,[0,0],'Color','black'); hold off;
%% Build structured T(fhat) matrix
addpath('../etc');
K0 = 3; %filter size (2*K0+1)x(2*K0+1) %trim Fourier data so that T is close to square
sampsiz = [11,11]; %necessary minimum number of samples
kx = repmat(-20:20,[41,1]);
ky = rot90(kx,1);
k(1,:) = kx(:); k(2,:) = ky(:);
sampinds = get_lowpass_inds(k,sampsiz);
fxhat_trim = reshape(fxhat(sampinds),sampsiz);
fyhat_trim = reshape(fyhat(sampinds),sampsiz);
Tx = im2col(rot90(fxhat_trim,2),[2*K0+1,2*K0+1]).';
Ty = im2col(rot90(fyhat_trim,2),[2*K0+1,2*K0+1]).';
%% Show necessary condition allows recovery of minimal polynomial
Tfull = [Tx(1:(end-1),:);Ty(1:(end-1),:)]; %remove last row of Tx and Ty to reach necessary limit
fprintf('Dimensions of T(fhat) = %d x %d\n',size(Tfull,1),size(Tfull,2));

% Compute SVD to obtain nullspace
[U,S,V] = svd(Tfull);
s = diag(S);
s = s/s(1);
logS1 = log10(s);
figure(3); stem(logS1); title('Normalized singular values of T(fhat) (log scale)')

% Plot filter cooresponding to minimal singular value
m=512;
ind = mod((1:(2*K0+1))-(K0+1),m)+1;
c1 = V(:,end);
F1 = zeros(m);
F1(ind,ind) = reshape(c1,[2*K0+1,2*K0+1]);
F1 = (m^2*((ifft2(F1))));
figure(4); 
imagesc(-real(F1)/max(real(F1(:))),[-0.5,0.5]); colormap parula; colorbar; title('recovery from necessary minimum of samples');
hold on;
contour(real(F1),[0 0],'Color','black');
hold off;
%% Show that for less than necessary condition the recovery fails
Tdef = [Tx(1:(end-2),:);Ty(1:(end-2),:)]; %remove last row two rows of Tx and Ty to go below necessary bound
fprintf('Dimensions of T(fhat) = %d x %d\n',size(Tdef,1),size(Tdef,2));

[U,S,V] = svd(Tdef);
s = diag(S);
s = s/s(1);
logS2 = log10(s);
figure(5); stem(logS2); title('Normalized singular values of rank deficient T(fhat) (log scale)')

% Plot filter cooresponding to minimal singular value
m=512;
ind = mod((1:(2*K0+1))-(K0+1),m)+1;
c2 = V(:,end);
F2 = zeros(m);
F2(ind,ind) = reshape(c2,[2*K0+1,2*K0+1]);
F2 = (m^2*((ifft2(F2))));
figure(6); 
imagesc(real(F2)/max(real(F2(:))),[-0.1,0.1]); colormap parula; colorbar; title('recovery from below necessary minimum samples');
hold on;
contour(real(F2),[0 0],'Color','black');
hold off;
%% 
% %% export plot: original image
% fig = figure(10);
% imagesc(fsmall,[-0.5,1]); colormap parula; colorbar;
% axis image
% axis off
% %export_fig edge_recovery_pwc_n3.eps -transparent
% %export_fig edge_recovery_pwc_n7.eps -transparent
% %% export plot: sucessful recovery/minimal polynomial
% fig = figure(11);
% scale = [-1, 1];
% mask = real(F1)/max(real(F1(:)));
% mask2 = sign(mask).*abs(mask).^(1/2);
% imagesc(mask,scale); colormap parula ; colorbar;
% axis image
% axis off
% hold on
% contour(real(F1),[0 0],'Color','black');
% hold off
% %export_fig edge_recovery_mask1_n3.eps -transparent
% %export_fig edge_recovery_mask1_n7.eps -transparent
% %% export plot: failed recovery
% fig = figure(12);
% mask = real(F2)/max(real(F2(:)));
% mask2 = sign(mask).*abs(mask).^(1/2);
% imagesc(-mask,[-0.05,0.05]); colormap parula ; colorbar;
% axis image
% axis off
% hold on
% contour(real(F2),[0 0],'Color','black');
% hold off
% %export_fig edge_recovery_mask2_n3.eps -transparent
% %export_fig edge_recovery_mask2_n7.eps -transparent
% %% plot: singular vals -- successful recovery
% fig = figure(13);
% stem(logS1,'.-');
% axis([1 50 -8 0])
% set(fig,'Position',[100 100 200 150]);
% %export_fig edge_recovery_singvals1_n3.eps -transparent
% %export_fig edge_recovery_singvals1_n7.eps -transparent
% %% plot: singular vals -- successful recovery
% fig = figure(14);
% stem(logS2,'.-');
% axis([1 50 -8 0])
% set(fig,'Position',[100 100 200 150]);
% %export_fig edge_recovery_singvals2_n3.eps -transparent
% %export_fig edge_recovery_singvals2_n7.eps -transparent
