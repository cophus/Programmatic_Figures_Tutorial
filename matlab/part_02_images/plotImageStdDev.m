function [] = plotImageStdDev(imageInput,intRange,fileNameBase)

% Colin Ophus - clophus@lbl.gov - 2020 February
% APS tutorial example - Save images with the intensity scaled by the
%                        standard deviation of the intensity distribution.
%                        If fileNameBase input is not provided, plot image.

% Inputs:
% imageInput - scalar image array
% fileNameBase - char string containing file name

% Convert to double floating point
imageInput = double(imageInput);

% default settings for intensity range
if nargin < 2 || isempty(intRange)
   intRange = [-2 2]; 
end

% scale image into units of standard deviation
imageInput(:) = imageInput - mean(imageInput(:));
imageInput(:) = imageInput / sqrt(mean(imageInput(:).^2));

% scale image output using desired intensity range
imageInput(:) = (imageInput - intRange(1)) / (intRange(2) - intRange(1));

% clamp output from 0 to 1
imageInput(:) = min(max(imageInput,0),1);

% Save or plot image
if nargin < 3
    figure(1)
    clf
    imagesc(imageInput)
    axis equal off
    colormap(inferno)
    set(gca,'position',[0 0 1 1])
    caxis([0 1])
else
    fileName = [fileNameBase '.png'];
    imwrite(round(imageInput*255)+1,inferno,fileName);
end



end
