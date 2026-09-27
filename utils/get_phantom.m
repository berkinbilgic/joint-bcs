function [ contrast ] = get_phantom( size_ )

% returns Shepp Logan phantoms with different contrast


N = [size_,size_];
im1 = phantom(N(1));
im = (im1<0) .* im1;
im1 = im1 - 2*im;

% find different gray levels in phantom
ind = find(im1 > 0, 1);
lvl = im1(ind);
lvls = lvl;

im2 = im1;
k = 1;

while 1
    IM(:,:,k) = (im2==lvl) .* im2;
    k = k + 1;
    im2 = im2 - (im2==lvl) .* im2;
    lvl = find(im2 > 0, 1) ;
    if isempty(lvl)
        break
    end
    lvl = im2(lvl);
    lvls(k) = lvl;
end



%lvls(3) = 5e-17

IM = IM>0;
if size_<=32

    contrast(:,:,2) = IM(:,:,1) * lvls(2) + IM(:,:,2) * lvls(5) + IM(:,:,3) * lvls(3)...
                  + IM(:,:,4) * lvls(1);
    contrast(:,:,3) = IM(:,:,1) * lvls(5) + IM(:,:,2) * lvls(2) + IM(:,:,3) * lvls(1)...
                  + IM(:,:,4) * lvls(3);
    contrast(:,:,4) = IM(:,:,1) * lvls(2) + IM(:,:,2) * lvls(1) + IM(:,:,3) * lvls(3)...
                  + IM(:,:,4) * lvls(4);
    contrast(:,:,1) = im1;
    show_ = make_tile(contrast,2,2);
    
else 
    contrast(:,:,2) = IM(:,:,1) * lvls(6) + IM(:,:,2) * lvls(5) + IM(:,:,3) * lvls(3)...
                  + IM(:,:,4) * lvls(2) + IM(:,:,5) * lvls(1) + IM(:,:,6) * lvls(4);
    contrast(:,:,3) = IM(:,:,1) * lvls(5) + IM(:,:,2) * lvls(2) + IM(:,:,3) * lvls(1)...
                  + IM(:,:,4) * lvls(3) + IM(:,:,5) * lvls(4) + IM(:,:,6) * lvls(6);
    contrast(:,:,4) = IM(:,:,1) * lvls(2) + IM(:,:,2) * lvls(1) + IM(:,:,3) * lvls(3)...
                  + IM(:,:,4) * lvls(6) + IM(:,:,5) * lvls(4) + IM(:,:,6) * lvls(5);
    contrast(:,:,5) = IM(:,:,1) * lvls(4) + IM(:,:,2) * lvls(2) + IM(:,:,3) * lvls(6)...
                  + IM(:,:,4) * lvls(5) + IM(:,:,5) * lvls(3) + IM(:,:,6) * lvls(1);
    contrast(:,:,6) = IM(:,:,1) * lvls(6) + IM(:,:,2) * lvls(4) + IM(:,:,3) * lvls(1)...
                  + IM(:,:,4) * lvls(5) + IM(:,:,5) * lvls(2) + IM(:,:,6) * lvls(3);              
    contrast(:,:,1) = im1;

    show_ = make_tile(contrast,2,3);
end

imsc(show_,1)

end

