function [ show_ ] = make_tile( imgs, row_num, col_num, im_norm )

% cat 2d images together as given in row_num, col_num
% if im_norm > 0 => normalize each image to [0, 1] before cat

if row_num * col_num ~= size(imgs,3)
    disp('In make_tile: sizes do not match')
end

show_ = zeros([size(imgs,1)*row_num, size(imgs,2)*col_num]);

if nargin < 4
    im_norm = 0;
end

for r = 1:row_num
    s = imgs(:,:,col_num*(r-1)+1);
    if im_norm>0
        s = abs(s);
        s = s - min(s(:));
        s = s / max(s(:));
    end
    
    S = s;
    
    
    for c = 2:col_num
        S = cat(2, S, imgs(:,:,col_num*(r-1)+c));
    end
    
    if r == 1
        show_ = S;
    else
        show_ = cat(1, show_, S);
    end
end

end

