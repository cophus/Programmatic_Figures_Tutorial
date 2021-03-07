function [Irgb] = colorComplex(amp,phase,ampRange)

% Make an RGB image from complex data with amplitude scaling
amp = (amp - ampRange(1))/(ampRange(2) - ampRange(1));
amp = min(max(amp,0),1);
N = size(amp);
Irgb = ones(N(1),N(2),3);
Irgb(:,:,1) = mod(-phase/(2*pi),1);
Irgb(:,:,3) = amp;
Irgb = hsv2rgb(Irgb);



end