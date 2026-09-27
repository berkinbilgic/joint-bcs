function [tmp] = prohibitSampling(pdf,prohibit,radius,tol)

% generate undersampling pattern that does not coincide with the samples in
% the variable 'prohibit'
% but fully sample inside a circle of radius r

%	pdf - probability density function to choose samples from
%	tol  - the deviation from the desired number of samples in samples
%   prohibit - do not put samples on these pixels

% returns:
%	tmp - sampling pattern


pdf(find(pdf>1)) = 1;
K = sum(pdf(:));


sx = size(prohibit,1);
sy = size(prohibit,2);

[x,y] = meshgrid(linspace(-1,1,sy),linspace(-1,1,sx));
r = sqrt(x.^2+y.^2);
r = r/max(abs(r(:)));			

idx = r<radius;
prohibit(idx) = 0;
%imshow(prohibit)
tmp = zeros(size(pdf));
tmp(idx) = 1;

while sum(tmp(:)) < K-tol
     tmp = tmp + (rand(size(pdf)) < .1*pdf);
     tmp(tmp>1) = 1;
     tmp(prohibit) = 0;
     disp(num2str(mean(tmp(:))))
%     imshow(tmp),pause(.1)
end

end


