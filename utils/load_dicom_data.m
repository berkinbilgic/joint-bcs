function [ contrast ] = load_dicom_data( new_size )

% returns clinical MRI scans, discarding the 3rd one (causes fft problem)

if nargin<1
    new_size = [256,256];
end

for n = 1:5
    if n ~= 3      % problem with img3.ima
        c = dicomread(['img',num2str(n),'.ima']);
    end
    
    if n == 2
        c = cat(2, zeros(256,8),c,zeros(256,8));
    end
    
    if n == 4
        c = cat(2, zeros(256,12),c,zeros(256,12));
    end
    
    
    c = imresize(c,new_size);
            
    if n > 3
        contrast(:,:,n-1) = double(c);
    elseif n < 3
        contrast(:,:,n) = double(c);
    end
end


contrast = contrast / max(contrast(:));


cshow = cat(2, contrast(:,:,1), contrast(:,:,2));
cshow = cat(1, cshow, cat(2,contrast(:,:,3), contrast(:,:,4)));
figure(1), imagesc(cshow.^1), axis image, colorbar, colormap(gray)


end

