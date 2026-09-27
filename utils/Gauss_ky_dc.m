function [ m2d ] = Gauss_ky_dc( size_data, perc, dc_sigma )

% try to generate sampling pattern along ky with a Gaussian pdf 
% fill the dc region with ones
% dc_sigma =>  make (dc_sigma/6) of the dc region ones

if nargin < 3
    dc_sigma = 1;
end


m2d = zeros(size_data);
dc_len =  dc_sigma * size_data(2) / 12 ;
mu = (1 + size_data(2))/2;

z = round(mu-dc_len:mu+dc_len)';
Sigma = size_data(2) / 5;  % make 5*sigma = image_width


while length(z) / size_data(2) < perc
    n_samp = round(size_data(2) * (perc/10));

    z = [z; round( mu + Sigma * randn(n_samp,1) )];
    
    z(z <= 0) = 1;
    z(z > size_data(2)) = size_data(2);
    
    z = unique(z);
end


m2d(:,z) = 1;
disp(['Undersampling ratio: ', num2str(mean(m2d(:)))])
imsc(m2d,1);

end

