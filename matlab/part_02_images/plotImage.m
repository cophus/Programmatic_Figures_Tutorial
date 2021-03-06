function [imageRGB] = plotImage(imageInput,varargin)


% Colin Ophus - clophus@lbl.gov - 2020 February
% APS tutorial example 01 - plotting images with automatic color scaling

% In this example, we will plot scalar (single-value per pixel) images,
% with automatically determined color ranges, using a defined colormap.

% Input:
% imageInput - the raw image data, 1 value per pixel
flagSimpleScaling = false;  % Set to true to use mean / standard deviation scaling
flagDrawColorbar = true;  % set to true to draw a colorbar

% Optional inputs
% 2 element vector [min_value max_value] specifies the relative or abs color range
% [N x 3] array will be assigned as the colormap
% 'absolute' character vector to change the color scaling to absolute 
% 'ordered' character vector to change the color scaling to order instead 
% 'asymmetric' character vector to change the color scaling to asymmetric 
% 'log' character vector to change the color scaling to absolute 
% 'zero' character vector to fix black level at zero 
% any other character vector will be assigned as file name for writing an output png file

% Init
intRange = [];
cmap = [];
flagAbsoluteScaling = false;
flagOrderedScaling = false;
flagAsymmetricScaling = false;
flagLogScaling = false;
flagZeroScaling = false;
flagWriteOutput = false;

% Parse inputs
numInputs = nargin - 1;
for a0 = 1:numInputs
    if ischar(varargin{a0}) == true
        
        
        % Check for keywords as inputs
        if strcmp(varargin{a0},'absolute')
            flagAbsoluteScaling = true;
        elseif strcmp(varargin{a0},'ordered')
            flagOrderedScaling = true;
        elseif strcmp(varargin{a0},'asymmetric')
            flagAsymmetricScaling = true;
        elseif strcmp(varargin{a0},'log')
            flagLogScaling = true;
        elseif strcmp(varargin{a0},'zero')
            flagZeroScaling = true;
        else
            % Otherwise assume string represents output file name
            fileNameBase = varargin{a0};
            flagWriteOutput = true;
        end
        
    else
        % If input is not a string, check for range or colormap
        if numel(varargin{a0}) == 2
            intRange = varargin{a0};
            
        elseif size(varargin{a0},1) > 1 && size(varargin{a0},2) == 3
            cmap = varargin{a0};
            
        else
            disp(['Input parameter #' num2str(a0+1) ' was ignored.']);
        end
    end
end

% Default parameters
if isempty(intRange)
    intRange = [-1 1]*2;
end
if isempty(cmap)
    cmap = inferno(256);
end

% Scale figure
imageScale = double(real(imageInput));
if flagAbsoluteScaling == true
    imageScale(:) = (imageScale - intRange(1)) / (intRange(2) - intRange(1));
    
elseif flagOrderedScaling == true
    numberColors = size(cmap,1);
    
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
    
    % Generate image
    imageRGB = zeros(imageSize(1),imageSize(2),3);
    for a0 = 1:(length(indsRange)-1)
        inds = indsOrder((indsRange(a0)+1):(indsRange(a0+1)));
        % red
        imageRGB(inds) = cmap(a0,1);
        % green
        imageRGB(inds+(prod(imageSize))) = cmap(a0,2);
        % blue
        imageRGB(inds+(2*prod(imageSize))) = cmap(a0,3);
    end  
    
elseif flagAsymmetricScaling == true
    imageMedian = median(imageScale(:));
    sub = imageScale < imageMedian;
    imageRMSneg = sqrt(mean((imageScale(sub) - imageMedian).^2));
    sub = imageScale > imageMedian;
    imageRMSpos = sqrt(mean((imageScale(sub) - imageMedian).^2));
    
    intRangeNew = imageMedian + ...
        [intRange(1)*imageRMSneg intRange(2)*imageRMSpos];
    imageScale(:) = (imageScale - intRangeNew(1)) ...
        / (intRangeNew(2) - intRangeNew(1));
    
elseif flagLogScaling == true
    % logarithmic scaling
    sub = imageScale > 0;
    sig = log(imageScale(sub));
    
    sigMean = mean(sig);
    sigRMS = sqrt(mean((sig - sigMean).^2));
    intRangeNew = exp(sigMean + intRange*sigRMS);
    
    imageScale(:) = (imageScale - intRangeNew(1)) ...
        / (intRangeNew(2) - intRangeNew(1));
    
elseif flagZeroScaling == true
    % Fix zero level, compute RMS intensity relative to zero
    imageScale = imageScale / sqrt(mean(imageScale(:).^2));
    imageScale(:) = imageScale / intRange(2);
    
else
    % Simple scaling by mean +/- standard deviation
    imageScale = imageScale - mean(imageScale(:));
    imageScale = imageScale / sqrt(mean(imageScale(:).^2));
    imageScale(:) = (imageScale - intRange(1)) / (intRange(2) - intRange(1));
    
end

if flagOrderedScaling == false
    % Clamp values to 0 and 1
    imageScale(:) = min(max(imageScale,0),1);
    
    % Generate RGB image
    imageRGB = ind2rgb(round(imageScale*(size(cmap,1)-1))+1,cmap);
end

% If needed, write figure into file.
% Descriptive file name
if flagWriteOutput == true
    if flagAbsoluteScaling == true
        fileName = [fileNameBase '_absolute_range_' ...
            sprintf('%.01d',intRange(1)) '_' ...
            sprintf('%.01d',intRange(2)) ...
            '.png'];
    elseif flagAsymmetricScaling == true
        fileName = [fileNameBase '_asymmetric_range_' ...
            sprintf('%.01d',intRange(1)) '_' ...
            sprintf('%.01d',intRange(2)) ...
            '.png'];
    elseif flagOrderedScaling == true
        fileName = [fileNameBase '_ordered.png'];
    elseif flagLogScaling == true
        fileName = [fileNameBase '_log_range_' ...
            sprintf('%.01d',intRange(1)) '_' ...
            sprintf('%.01d',intRange(2)) ...
            '.png'];
    elseif flagZeroScaling == true
        fileName = [fileNameBase '_zero_range_' ...
            sprintf('%.01d',intRange(1)) '_' ...
            sprintf('%.01d',intRange(2)) ...
            '.png'];
    else
        fileName = [fileNameBase '_range_' ...
            sprintf('%.01d',intRange(1)) '_' ...
            sprintf('%.01d',intRange(2)) ...
            '.png'];
    end
    
    imwrite(imageRGB,fileName,'png')
end


% Plot the figure will full axes
figure(1)
clf
imagesc(imageRGB)
axis equal off
set(gca,'position',[0 0 1 1])



if flagDrawColorbar == true && flagOrderedScaling == true
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
end


end