function [  ] = tile( imgs, row_num, col_num, fig_num, color_bar, im_norm )

% cat 2d images together as given in row_num, col_num
% if im_norm > 0 => normalize each image to [0, 1] before cat

if row_num * col_num ~= size(imgs,3)
    disp('In make_tile: sizes do not match')
end

show_ = zeros([size(imgs,1)*row_num, size(imgs,2)*col_num]);

if nargin < 6
    im_norm = 0;
end
if nargin < 5
    color_bar = 0;
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

if nargin < 4
    imsc(show_)
else
    imsc(show_,fig_num,color_bar)
end

end

