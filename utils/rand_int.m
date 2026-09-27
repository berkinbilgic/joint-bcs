function [ out ] = rand_int( maxi, N )

% returns N unique random integers between 1 and maxi

if N > maxi
    disp('Error in rand_int')
    out = 0;
else
    step = round(N/10);
    out = unique( unidrnd(maxi,N,1) );
    while length(out) < N
        out = unique([out; unidrnd(maxi,step,1)]);
    end
    
    if length(out)>N
        out(N+1:end) = [];
    end
end


end

