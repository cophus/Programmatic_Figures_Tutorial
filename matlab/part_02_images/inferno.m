function [cmap] = inferno(numberColors)

% Colin Ophus - clophus@lbl.gov - 2020 February
% Quick implementation / tweak of a perceptually-uniform colormap "inferno" 
% by Stefan van der Walt and Nathaniel Smith.  If you run this script
% without an output, it will plot the colormap in a figure window.

% Input:
% numberColors - how many colors to include (# of RGB values). Default: 256

% Output:
% [numberColors x 3] sized array containing the colormap, columns = [R G B]

% default number of colors
if nargin < 1
    numberColors = 256;
end

% Init
cmap = ones(numberColors,3);
cmap(:,3) = linspace(0,1,numberColors);

% Determine specific color indices
indsEdge = round(numberColors*[0 0.3 0.5 0.75 0.90 1]) + [1 0 0 0 0 0];

% black -> purple
inds = indsEdge(1):(indsEdge(2)-1);
c = linspace(0,1,indsEdge(2)-indsEdge(1)+1)'; 
c(end) = [];
cmap(inds,1) = 4/6 + c/6;
cmap(inds,2) = 3*c.^2 - 2*c.^3;

% purple -> red
inds = indsEdge(2):(indsEdge(3)-1);
c = linspace(0,1,indsEdge(3)-indsEdge(2)+1)'; 
c(end) = [];
cmap(inds,1) = 5/6 + c/6;

% red -> orange
inds = indsEdge(3):(indsEdge(4)-1);
c = linspace(0,1,indsEdge(4)-indsEdge(3)+1)'; 
c(end) = [];
cmap(inds,1) = c/12;

% orange -> yellow
inds = indsEdge(4):(indsEdge(5)-1);
c = linspace(0,1,indsEdge(5)-indsEdge(4)+1)'; 
c(end) = [];
cmap(inds,1) = 1/12 + c/12;

% yellow -> white
inds = indsEdge(5):indsEdge(6);
c = linspace(0,1,indsEdge(6)-indsEdge(5)+1)'; 
cmap(inds,1) = 1/6;
cmap(inds,2) = 1 - 3*c.^2 + 2*c.^3;


% convert from HSV to RGB
cmap(:) = hsv2rgb(cmap);

% plotting
if nargout < 1
    figure(101)
    clf
    imagesc(reshape(cmap,[numberColors 1 3]))
    axis off
    set(gca,'xtick',[])

end