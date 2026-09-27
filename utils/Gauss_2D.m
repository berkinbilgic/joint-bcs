function [ m2d ] = Gauss_2D( size_data, perc )

% tries to implement sampling in 2d with Gaussian pdf

len = size_data(1);
for n = 2:length(size_data)
    len = len * size_data(n);
end

k = 1;
mu = (1 + size_data)/2;
Sigma = [5*size_data(1) 0 ; 0 5*size_data(2)];  % this is too emprical...
R = chol(Sigma);

while 1
    n_samp = round(len * (perc/5));

    if k == 1
        k = 2;
        z = round( repmat(mu,n_samp,1) + randn(n_samp,2)*R );
    else
        z = [z; round( repmat(mu,n_samp,1) + randn(n_samp,2)*R )];
    end
    
    for n = 1:2
        z_ = z(:,n);
        z_(z_ <= 0) = 1;
        z_(z_ > size_data(n)) = size_data(n);
        z(:,n) = z_;
    end

    z = unique(z,'rows');    
    
    if size(z,1) / len >= perc
        break
    end
end

disp(['Undersampling ratio: ', num2str(size(z,1) / len)])

m2d = zeros(size_data);
for n = 1:size(z,1)    
    m2d(z(n,1),z(n,2)) = 1;
end


end