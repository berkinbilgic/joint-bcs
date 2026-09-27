function [ wav ] = db4( img, num_scale )

% try to get db4 wavelet at 'num_scale' scales

siz_ = log2(size(img));

if sum(siz_ - round(siz_)) ~= 0
    disp('In db4: input size not a power of 2')
    return
end


[Lo,Hi] = compute_wavelet_filter('Daubechies',4);  % Db4


wav = img;
wav_sub = wav;
for n = 1:num_scale

    Coarse = subsampling( cconv(wav_sub,Lo,1),1 );
    Detail = subsampling( cconv(wav_sub,Hi,1),1 );
    wav_sub = cat(1, Coarse, Detail );


    Coarse = subsampling( cconv(wav_sub,Lo,2),2 );
    Detail = subsampling( cconv(wav_sub,Hi,2),2 );
    wav_sub = cat(2, Coarse, Detail );

    wav(1:size(wav_sub,1), 1:size(wav_sub,2)) = wav_sub;
    
    wav_sub = wav(1:size(wav_sub,1)/2, 1:size(wav_sub,2)/2);
    
end


end

