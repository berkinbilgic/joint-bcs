function [ A,B ] = get_dft_matrix( size_ )

% returns dft matrices along x (A) and y (B)

N = size_(1);
M = size_(2);

A = sparse([], [], [], N*M, N*M, (N^2)*M);

Mat = zeros(N);
wN = exp(-2*pi*1i/N);
for n = 1:N
    Mat(n,:) = wN.^((0:N-1)*(n-1));
end

for n = 1:M
    A(N*(n-1)+1:N*n, N*(n-1)+1:N*n) = Mat;
end


B = sparse([], [], [], N*M, N*M, (M^2)*N);
wM = exp(-2*pi*1i/M);

wMat = wM.^(0:M-1);
disp('Generating dft_y matrix: ')

for k = 1:M
    temp = wMat.^(k-1);
    disp([num2str(100*k/M),' % complete'])
    for n = 1:N
        B( n + (k-1) * N, n:N:N*M ) = temp;
    end
end


end

