function [imageRGB] = colorComplex(EW,ampRange,ampPower)

% Colin Ophus - 2021 March
% HSV version of exit wave coloring

imageRGB = ones(size(EW,1),size(EW,2),3);
imageRGB(:,:,1) = mod(angle(EW)/(2*pi),1);
imageRGB(:,:,3) = min(max( ....
     (abs(EW) - ampRange(1)) / (ampRange(2) - ampRange(1)),0),1);
 if nargin > 2
    imageRGB(:,:,3) = imageRGB(:,:,3).^ampPower;
end
imageRGB(:) = hsv2rgb(imageRGB);

if nargout == 0
    figure(11)
    clf
    imagesc(imageRGB)
    axis equal off
    set(gca,'position',[0 0 1 1])
end

end