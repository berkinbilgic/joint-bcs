function [ TVx, TVs ] = get_tv( IXFMx, IXFMs )

Dx = IXFMx([2:end,end],:,:) - IXFMx;
Dy = IXFMx(:,[2:end,end],:) - IXFMx; 
TVx = cat(3,Dx,Dy);

TVs = 0;
if IXFMs ~= 0
    Dx = IXFMs([2:end,end],:,:) - IXFMs;
    Dy = IXFMs(:,[2:end,end],:) - IXFMs; 
    TVs = cat(3,Dx,Dy);
end

end

