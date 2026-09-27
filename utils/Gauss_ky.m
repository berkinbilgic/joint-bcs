function [ m2d ] = Gauss_ky( size_data, perc )

% try to generate sampling pattern along ky with a Gaussian pdf 

m2d = zeros(size_data);


k = 1;
mu = (1 + size_data(2))/2;
Sigma = size_data(2) / 6;  % make 6*sigma = image_width


while 1
    n_samp = round(size_data(2) * (perc/8));

    if k == 1
        k = 2;
        z = round( mu + Sigma * randn(n_samp,1)  );
    else
        z = [z; round( mu + Sigma * randn(n_samp,1) )];
    end
    
    z(z <= 0) = 1;
    z(z > size_data(2)) = size_data(2);
    
    z = unique(z);
    
    if length(z) / size_data(2) >= perc
        break
    end
end

disp(['Undersampling ratio: ', num2str(length(z) / size_data(2))])
m2d(:,z) = 1;
imsc(m2d,1);

end

