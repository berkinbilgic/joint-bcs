function [ wav ] = inv_db4( wav, num_scale )

% try to get inverse db4 wavelet at 'num_scale' scales


[Lo,Hi] = compute_wavelet_filter('Daubechies',4);  % Db4

start_size = size(wav) / 2^(num_scale-1);

for n = 1:num_scale

    w = wav(1:start_size(1),1:start_size(2));
    Coarse = w(1:size(w)/2,:);
    Detail = w(1+size(w)/2:end,:);
    
    Coarse = cconv( upsampling(Coarse,1), reverse(Lo),1 );
    Detail = cconv( upsampling(Detail,1), reverse(Hi),1 );
    A = Coarse + Detail;

    Coarse = A(:,1:size(A,2)/2);
    Detail = A(:,1+size(A,2)/2:end);
    
    Coarse = cconv( upsampling(Coarse,2), reverse(Lo),2 );
    Detail = cconv( upsampling(Detail,2), reverse(Hi),2 );
    A = Coarse + Detail;
   
    wav(1:size(w,1),1:size(w,2)) = A;

    start_size = start_size * 2;
end


end

