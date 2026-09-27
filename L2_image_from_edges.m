function X = L2_image_from_edges(Xhatp,DvX,DhX,beta)

%This function will attempt to recover an image (hopefully sparse in
%gradient) from its edge estimates by solving a least squares problem of
%the form

%X = argmin_Z ||Zv - DvX||_2^2 + ||Zh - DhZ||_2^2 + lambda*||Zhat|_{Omega}
%- Xhatp||_2^2

%where lambda will be zero outside of Omega and effectively infinity on Omega.


N = size(Xhatp,1);
W = zeros(N,N);
W(find(Xhatp)) = 1;

[k2,k1] = meshgrid(0:N-1,0:N-1);

FDV = (1 - exp(-2*pi*i*k1/N));
FDH = (1 - exp(-2*pi*i*k2/N));
EH = zeros(N,N);
EH(1,1) = 1;

%Xhat = ((conj(FDV).*fft2(DvX) + conj(FDH).*fft2(DhX))./((1+beta)*(abs(FDV).^2+abs(FDH).^2+EH)));
%Xhat(1,1) = Xhatp(1,1);
%X = ifft2(Xhat);

X = ifft2(((conj(FDV).*fft2(DvX) + conj(FDH).*fft2(DhX))./((1+beta)*(abs(FDV).^2+abs(FDH).^2+EH))).*(1-W) + Xhatp);
