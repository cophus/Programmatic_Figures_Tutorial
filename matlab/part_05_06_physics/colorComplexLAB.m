function [imageRGB] = colorComplexLAB(EW,ampRange,ampPower)

% Colin Ophus - 2021 March
% L*a*b version of exit wave coloring

% abRange = 64;

p = angle(EW);
amp = abs(EW);
amp = min(max((amp - ampRange(1)) / (ampRange(2) - ampRange(1)),0),1);
if nargin > 2
    amp = amp.^ampPower;
end

imageRGB = ones(size(EW,1),size(EW,2),3);
imageRGB(:,:,1) = mod(p/(2*pi),1);
imageRGB(:,:,3) = amp;
imageRGB(:) = hsv2rgb(imageRGB);


imageRGB(:) = rgb2lab(imageRGB);
% min(min(imageRGB(:,:,1)))
imageRGB(:,:,1) = 100 * amp;
imageRGB(:) = lab2rgb(imageRGB);

% % rgb = lab2rgb(lab)
% % abRange = 200;
% % t = repmat(linspace(-abRange,abRange,size(imageRGB,2)),[size(imageRGB,1) 1]);
% % imageRGB(:,:,2) = t;
% imageRGB = zeros(size(EW,1),size(EW,2),3);
% imageRGB(:,:,1) = 100 * amp;
% imageRGB(:,:,2) = cos(p) * abRange;
% imageRGB(:,:,3) = sin(p) * abRange;
% imageRGB(:) = lab2rgb(imageRGB);


if nargout == 0
    figure(11)
    clf    
    imagesc(imageRGB)
    axis equal off
    set(gca,'position',[0 0 1 1])
end


end
