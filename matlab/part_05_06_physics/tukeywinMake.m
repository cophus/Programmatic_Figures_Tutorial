function [windowOutput] = tukeywinMake(windowSize,windowFrac)

% Colin Ophus - 2020 Sept
% Quick re-definition of Tukey / Hann window function


if nargin == 1
    windowFrac = [0 0];
elseif length(windowSize) == 2 && length(windowFrac) == 1
    windowFrac = [1 1]*windowFrac;
end

% 1D window
windowOutput = min(1/(1 - windowFrac(1)) * (1 - ...
    abs((windowSize(1)+1)/2 - (1:windowSize(1))') * 2 / windowSize(1)),1);
windowOutput(:) = sin(windowOutput*(pi/2)).^2;

% 2D window
if length(windowSize) == 2
    wy = min(1/(1 - windowFrac(2)) * (1 - ...
        abs((windowSize(2)+1)/2 - (1:windowSize(2))) * 2 / windowSize(2)),1);
    wy(:) = sin(wy*(pi/2)).^2;
    
    windowOutput = windowOutput * wy;
end


end