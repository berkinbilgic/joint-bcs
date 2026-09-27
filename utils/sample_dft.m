function [ dft, m2d_orig, ind ] = sample_dft( size_data, pdf, perc, dft )

%function [ dft, m2d, ind, m2d_orig ] = sample_dft( size_data, pdf, perc, dft )

% try to map shifted fft sampling mask to not-shifted dft matrix
% perc = ratio of k-space samples we want to keep

% Gaussian mask

if strcmp(pdf,'Gaussian')
    m2d = Gauss_2D( size_data, perc );
else
    if strcmp(pdf,'Uniform') || strcmp(pdf,'uniform')
        m2d = rand(size_data) > (1-perc);
    else
        disp('Error in sample_dft')
    end
end

m2d_orig = m2d;
m2d = fftshift(m2d);
m2d = 1-m2d;
ind = find(m2d);

dft(ind,:) = [];

end

