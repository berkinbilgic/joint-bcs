function [ x ] = CG( x0, mask, y, param )

switch param.xfm_type
    case 0
        % only Identity transform
        param.tv_weight = 0;
     case 1
        % only TV norm
        param.xfm_weight = 0;
    case 2
        % only wavelet transform
        param.tv_weight = 0;
    case 3
        % only curvelet transform, using default parameters
        param.tv_weight = 0;
    case 4
        % identity & TV norm
    case 5
        % wavelet & TV norm
    case 6
        % curvelet & TV norm
end


fft_coef = sqrt(length(mask(:)));

IXFMx = get_invXFM(x0, 0, param);
X = get_fft(IXFMx, 0, mask, fft_coef);
TVx = get_tv(IXFMx, 0, param);


g_o = grad_CG( X, x0, TVx, y, fft_coef, param );
s = -g_o;
x = x0;
g_new = g_o;
dot_new = dot(g_new(:),g_new(:));
t0 = 1;


for k = 1:param.max_iter

    % line search
    alfa_ = t0;
    
    % precompute fft's 
    [IXFMx, IXFMs] = get_invXFM(x, s, param);
    [X,S] = get_fft(IXFMx, IXFMs, mask, fft_coef);
    [TVx,TVs] = get_tv(IXFMx, IXFMs,param);
    
    func_LHS = cost_function_CG(X + alfa_*S, x + alfa_*s, TVx + alfa_*TVs, y, param);
    func_RHS1 = cost_function_CG(X, x, TVx, y, param);
    func_RHS2 = .01 * abs( dot(g_new(:),s(:)) );
    
    LS_iter = 0;
    
    while func_LHS > func_RHS1 + alfa_ * func_RHS2  && LS_iter < 50
        alfa_ = .6 * alfa_;
        func_LHS = cost_function_CG(X + alfa_*S, x + alfa_*s, TVx + alfa_*TVs, y, param);
        LS_iter = LS_iter + 1;
    end
    
    
    if LS_iter > 2
        t0 = t0 * .6;
    elseif LS_iter < 1
		t0 = t0 / .6;
	end
    

    dot_o = dot_new;
    
    x = x + alfa_ * s;
    g_new = grad_CG( X + alfa_*S, x, TVx + alfa_*TVs, y, fft_coef, param );
    
    dot_new = dot(g_new(:),g_new(:));
    
    beta_ = dot_new / dot_o;
    s = -g_new + beta_ * s; 


    if( norm(s(:)) < 1e-2 )
        break
    end
    
end



end

% -------------------------------------------------------------------------
function [ g ] = grad_CG( X, x, TVx, y, fft_coef, param )

mu_ = 1e-15;

g = 2 * ifftshift( ifft2( fftshift( X - y ) ) ) * fft_coef;
if param.tv_weight > 0
    g = g + param.tv_weight * inv_tv( TVx.*(TVx.*conj(TVx) + mu_).^(-.5) );
end

if param.xfm_type == 2 || param.xfm_type == 5
   % wavelet
   g = g + param.xfm_weight * inv_db4( x .* ((x.*conj(x) + mu_).^(-.5)), param.num_scales);
else
    if param.xfm_type == 3 || param.xfm_type == 6
       % curvelet
       g = g + param.xfm_weight * ifdct_wrapping( x .* (x.*conj(x) + mu_).^(-.5) );  % need to take care of cell structure
    else
       % identity
       g = g + param.xfm_weight * x .* (x.*conj(x) + mu_).^(-.5);
    end
end


end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
function [ cost ] = cost_function_CG( X, x, TVx, y, param )

mu_ = 1e-15;

temp = X - y;
cost = dot(temp(:), temp(:));

xfm = 0;
if param.xfm_weight > 0
    xfm = sum( (x(:).*conj(x(:))+mu_).^.5 ) * param.xfm_weight; 
end

TV = 0;
if param.tv_weight > 0
    TV = sum( (TVx(:).*conj(TVx(:))+mu_).^.5 ) * param.tv_weight; 
end

cost = cost + TV + xfm;

end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
function [ IXFMx, IXFMs ] = get_invXFM( x, s, param )
% returns the inverse transform (converts to image space)

    IXFMs = 0;
    if param.xfm_type == 2 || param.xfm_type == 5
        % wavelet
        IXFMx = inv_db4(x,param.num_scales);
        if s ~= 0
            IXFMs = inv_db4(s,param.num_scales);
        end
    else
        if param.xfm_type == 3 || param.xfm_type == 6
           % curvelet
           IXFMx = ifdct_wrapping(x);
           if s ~= 0
              IXFMs = ifdct_wrapping(s);
           end
        else
            % identity transform
            IXFMx = x;
            if s ~= 0
               IXFMs = s;
            end 
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
function [ TVx, TVs ] = get_tv( IXFMx, IXFMs, param )
  if param.tv_weight == 0
      TVx = 0;
      TVs = 0;
  else
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