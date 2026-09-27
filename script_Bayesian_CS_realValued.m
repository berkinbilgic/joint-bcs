%% load multicontrast phantoms

addpath utils\

load new_SL

im_size = [128,128];
L = size(img,3);

tile(img,1,L,1)


%% undersample in 1-dimension

R = 2;
pdf = genPDF([im_size(2),1], 5, 1/R, 2, 0.08, 0);
m2d = [];

for t = 1:L
    m1d = fftshift(genSampling(pdf,10,1));
    m2d(:,:,t) = repmat(m1d,[im_size(1),1]);     
end

tile(fftshift(m2d),1,L,1)


%% undersample with power law pdf in 2-dimensions

R = 9;
pdf = genPDF([im_size,im_size], 5, 1/R, 2, 0.1, 0);	

for t = 1:L
    m2d(:,:,t) = fftshift(genSampling(pdf,10,60)); 
end
tile(fftshift(m2d),1,L,1)


%% compute k-space of image gradients

[k2,k1] = meshgrid(0:im_size(2)-1,0:im_size(1)-1);
fdx = 1 - exp(-2*pi*1i*k1/im_size(1));
fdy = 1 - exp(-2*pi*1i*k2/im_size(2));

fft_scale = sqrt(prod(im_size));


y = [];     zf = [];        loc = [];       tx = [];        ty = [];    img_x = [];     img_y = [];
for h = 1:L
    y(:,:,h) = fft2(img(:,:,h)).*m2d(:,:,h) / fft_scale;
    zf(:,:,h) = ifft2(y(:,:,h)) * fft_scale;
    
    loc{h} = find(y(:,:,h));
    y_dx = y(:,:,h).*fdx;
    y_dy = y(:,:,h).*fdy;
    
    tx{h} = cat(1, real(y_dx(loc{h})), imag(y_dx(loc{h})) );
    ty{h} = cat(1, real(y_dy(loc{h})), imag(y_dy(loc{h})) );

    img_x(:,:,h) = ifft2(fft2(img(:,:,h)).*fdx);
    img_y(:,:,h) = ifft2(fft2(img(:,:,h)).*fdy);    
end

tile(real(zf),1,L,1), title(['Zero filling RMSE: ', num2str(100*norm(zf(:)-img(:))/norm(img(:))), ' percent'])


%% joint reconstruction

max_iters = 1e4;

tic, dx = mt_CSfft2_rect(im_size, tx, loc, 0, 0, 1e-8, img_x, max_iters);      

     dy = mt_CSfft2_rect(im_size, ty, loc, 0, 0, 1e-8, img_y, max_iters);      toc

   

% Least Squares reconstruction from the gradient estimates

P1x = reshape(dx,[im_size,L]);
P1y = reshape(dy,[im_size,L]);

Phat = y * fft_scale;

Pls = [];
for h = 1:L
    Pls(:,:,h) = L2_image_from_edges_rect(Phat(:,:,h),P1x(:,:,h),P1y(:,:,h),0);
end

tile(img,1,L,1)
tile(real(Pls),1,L,2), title(['Bayesian CS RMSE: ', num2str(100*norm((Pls(:))-img(:))/norm(img(:))), ' percent'])


