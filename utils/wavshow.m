function [  ] = wavshow( wav, figno )

if nargin < 2
    imsc( abs(wav).^.4 )
else
    imsc( abs(wav).^.4 ,figno)
end

end

