function [v,basis_v] = basis_complex_rect(k, loc, M, im_size, fft_scale)
% try to generate the k^th basis vector in the dft matrix with size NxM
% where 'loc' contains the indices of the vector where the rows are not
% removed (due to undersampling)


%im_size = sqrt(M);

x = zeros(im_size);
x(k) = 1;

%basis_v = fft2(x) / sqrt(M);
basis_v = fft2(x) / fft_scale;

basis_v = basis_v(:);
basis_v(setdiff(1:M,loc)) = [];

v = [real(basis_v);imag(basis_v)];
end