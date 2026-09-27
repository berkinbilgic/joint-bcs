function [  ] = imsc( img, fig_num, colorbar_on, title_ )

if nargin<3
    colorbar_on = 1;
end

if nargin<2
    figure
else
    figure(fig_num)
end

if sum(abs(imag(img(:)))) > 0 
    img = abs(img);
end

if colorbar_on > 0
    imagesc(img), colorbar, colormap(gray), axis image
else
    imagesc(img), colormap(gray), axis image
end


if nargin == 4
    title(title_)
end

end

