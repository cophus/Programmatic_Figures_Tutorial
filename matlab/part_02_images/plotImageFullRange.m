function [] = plotImageFullRange(imageInput,fileNameBase)

% Colin Ophus - clophus@lbl.gov - 2020 February
% APS tutorial example - Save images from min to max range, inferno colormap
%                        Note that if fileNameBase is not specified, this
%                        script will plot the image instead, and no file 
%                        will be written.

% Inputs:
% imageInput - scalar image array
% fileNameBase - char string containing file name


% Convert to double floating point
imageInput = double(imageInput);

% scale image from 0 to 1
imageInput(:) = imageInput - min(imageInput(:));
imageInput(:) = imageInput / max (imageInput(:));

% Save or plot image
if nargin < 2
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
