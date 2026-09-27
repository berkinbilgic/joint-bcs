function [ x ] = CG_TV_XFM( x0, mask, y, param )

switch param.xfm_type
    case 0
        % only Identity transform
        tv_weight = 0;
        xfm_weight = param.xfm_weight;
        levels = 0;
    case 1
        % only TV norm
        tv_weight = param.tv_weight;
        xfm_weight = 0;
        levels = 0;
    case 2
        % only wavelet transform
        tv_weight = 0;
        xfm_weight = param.xfm_weight;
        levels = param.num_scales;    % no of wavelet scales
    case 3
        % only curvelet transform, using default parameters
        tv_weight = 0;
        xfm_weight = param.xfm_weight;
    case 4
        % identity & TV norm
        tv_weight = param.tv_weight;
        xfm_weight = param.xfm_weight;
        levels = 0;
    case 5
        % wavelet & TV norm
        tv_weight = param.tv_weight;
        xfm_weight = param.xfm_weight;
        levels = param.num_scales; 
    case 6
        % curvelet & TV norm
        tv_weight = param.tv_weight;
        xfm_weight = param.xfm_weight;
end


fft_coef = sqrt(length(x0(:)));

IXFMx = get_invwav(x0, 0, xfm_weight, levels);
X = get_fft(IXFMx, 0, mask, fft_coef);
TVx = get_tv(IXFMx, 0);

g_o = grad_CG( X, x0, TVx, y, xfm_weight, tv_weight, fft_coef, levels );
s = -g_o;
x = x0;
g_new = g_o;
dot_new = dot(g_new(:),g_new(:));
t0 = 1;


for k = 1:param.max_iter

    % line search
    alfa_ = t0;
    
    % precompute fft's 
    [IXFMx, IXFMs] = get_invwav(x, s, xfm_weight, levels);
    [X,S] = get_fft(IXFMx, IXFMs, mask, fft_coef);
    [TVx,TVs] = get_tv(IXFMx, IXFMs);
    
    func_LHS = cost_function_CG(X + alfa_*S, x + alfa_*s, TVx + alfa_*TVs, y, xfm_weight, tv_weight);
    func_RHS1 = cost_function_CG(X, x, TVx, y, xfm_weight, tv_weight);
    func_RHS2 = .01 * abs( dot(g_new(:),s(:)) );
    
    LS_iter = 0;
    
    while func_LHS > func_RHS1 + alfa_ * func_RHS2  && LS_iter < 50
        alfa_ = .6 * alfa_;
        func_LHS = cost_function_CG(X + alfa_*S, x + alfa_*s, TVx + alfa_*TVs, y, xfm_weight, tv_weight);
        LS_iter = LS_iter + 1;
    end
    
    
    if LS_iter > 2
        t0 = t0 * .6;
    elseif LS_iter < 1
		t0 = t0 / .6;
	end
    

    dot_o = dot_new;
    
    x = x + alfa_ * s;
    g_new = grad_CG( X + alfa_*S, x, TVx + alfa_*TVs, y, xfm_weight, tv_weight, fft_coef, levels );
    
    dot_new = dot(g_new(:),g_new(:));
    
    beta_ = dot_new / dot_o;
    s = -g_new + beta_ * s; 


    if( norm(s(:)) < 1e-2 )
        break
    end
    
end



end

% -------------------------------------------------------------------------
function [ g ] = grad_CG( X, x, TVx, y, XFM_weight, TV_weight, fft_coef, levels )

mu_ = 1e-15;

if XFM_weight > 0
    g = 2 * db4( ifftshift(ifft2(fftshift( X - y ))) * fft_coef, levels ); 
    g = g + XFM_weight * x .* (x.*conj(x) + mu_).^(-.5);
    g = g + TV_weight * db4( inv_tv( TVx.*(TVx.*conj(TVx) + mu_).^(-.5) ), levels );
else
    g = 2 * ifftshift(ifft2(fftshift( X - y))) * fft_coef; 
    g = g + TV_weight * inv_tv( TVx.*(TVx.*conj(TVx) + mu_).^(-.5) );
end

end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
function [ cost ] = cost_function_CG( X, x, TVx, y, XFM_weight, tv_weight )

mu_ = 1e-15;

temp = X - y;
cost = dot(temp(:), temp(:));

xfm = 0;
if XFM_weight > 0
    xfm = sum( (x(:).*conj(x(:))+mu_).^.5 ) * XFM_weight; 
end

TV = (TVx(:).*conj(TVx(:))+mu_).^.5; 
cost = cost + tv_weight * sum(TV(:)) + xfm;

end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
function [ IXFMx, IXFMs ] = get_invwav( x, s, XFM_weight, levels )

    if XFM_weight > 0
        IXFMx = inv_db4(x, levels);
    else
        IXFMx = x;
    end

    IXFMs = 0;
    if s ~= 0
        if XFM_weight > 0
            IXFMs = inv_db4(s, levels);
        else
            IXFMs = s;
        end
    end

end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
function [ X, S ] = get_fft( IXFMx, IXFMs, mask, fft_coef )

    X = fftshift(fft2(ifftshift(IXFMx))) .* mask / fft_coef;
    
    S = 0;
    if IXFMs ~= 0
        S = fftshift(fft2(ifftshift(IXFMs))) .* mask / fft_coef;
    end
end

% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
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
  
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
function res = inv_tv(tv_x)

re = tv_x([1,1:end-1],:,1) - tv_x(:,:,1);
re(1,:) = -tv_x(1,:,1);
re(end,:) = tv_x(end-1,:,1);

res = tv_x(:,[1,1:end-1],2) - tv_x(:,:,2);
res(:,1) = -tv_x(:,1,2);
res(:,end) = -tv_x(:,end-1,2);

res = res + re;
end
% -------------------------------------------------------------------------