function [imageRGB,intensityRange] = plotImageOrdered(imageInput,fileNameBase)

% Colin Ophus - clophus@lbl.gov - 2020 February
% APS tutorial example - Save images with the intensity scaled by the
%                        ordered distribution of intensties

% Inputs:
% imageInput - scalar image array
% fileNameBase - char string containing file name
flagDrawColorbar = true;  % draw and label a colorbar for the output image

% Output:
% imageRGB - output image with colormap applied

% Colormap
cmap = inferno;
% cmap = violetFire;
numberColors = size(cmap,1);

% Convert to double floating point
imageInput = double(imageInput);

% Ordering of pixel intensities
[sigOrder,indsOrder] = sort(imageInput(:));

% indices
imageSize = size(imageInput);
indsRange = round(linspace(0,1,numberColors+1)' * prod(imageSize));
if nargout > 1 || flagDrawColorbar == true
    % Output the intensity ranges corresponding to each color
    intensityRange = [ ...
        sigOrder(max(indsRange(1:end-1),1)) ...
        sigOrder(indsRange(2:end))];
end

% init
imageRGB = zeros(imageSize(1),imageSize(2),3);

% Generate image
for a0 = 1:(length(indsRange)-1)
    inds = indsOrder((indsRange(a0)+1):(indsRange(a0+1)));
    
    % red
    imageRGB(inds) = cmap(a0,1);
    
    % green
    imageRGB(inds+(prod(imageSize))) = cmap(a0,2);
    
    % blue
    imageRGB(inds+(2*prod(imageSize))) = cmap(a0,3);
end


% Save or plot image
if nargin < 2
    figure(1)
    clf
    imagesc(imageRGB)
    axis equal off
    set(gca,'position',[0 0 1 1])
else
    fileName = [fileNameBase '.png'];
    imwrite(imageRGB,fileName);
end


if flagDrawColorbar == true
    figure(2)
    clf
    set(gcf,'color','w')
    
    % Vertical orientation
    imagesc(reshape(cmap,[numberColors 1 3]));
    set(gca,'position',[0.1 0.05 0.1 0.90])
    set(gca,'xtick',[])
    set(gca,'ydir','normal')
    set(gca,'yaxislocation','right')
    % labels
    numTicks = 16 + 1;
    indsPlot = round(linspace(1,numberColors,numTicks));
    indValues = (intensityRange(indsPlot,1) + intensityRange(indsPlot,1))/2;
    indValues(:) = round(indValues);
    set(gca,'ytick',indsPlot);
    set(gca,'yticklabels',indValues);
    
    %     % horizontal orientation
    %     imagesc(reshape(cmap,[1 numberColors 3]));
    %     set(gca,'position',[0.05 0.45 0.90 0.10])
    %     set(gca,'ytick',[])
    %     % labels
    %     numTicks = 8 + 1;
    %     indsPlot = round(linspace(1,numberColors,numTicks));
    %     set(gca,'xtick',indsPlot);
    %     set(gca,'xticklabels', ....
    %         (intensityRange(indsPlot,1) + intensityRange(indsPlot,1))/2);
    
    
end

end