function res = inv_tv(tv_x, W)

tv_x = tv_x.*W;

re = tv_x([1,1:end-1],:,1) - tv_x(:,:,1);
re(1,:) = -tv_x(1,:,1);
re(end,:) = tv_x(end-1,:,1);

res = tv_x(:,[1,1:end-1],2) - tv_x(:,:,2);
res(:,1) = -tv_x(:,1,2);
res(:,end) = -tv_x(:,end-1,2);

res = res + re;

end



