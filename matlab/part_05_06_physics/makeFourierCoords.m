function [qx,qy] = makeFourierCoords(N,pixelSize)

% Colin Ophus - clophus@lbl.gov - 2020 February
% Generate Fourier space coordinates in x and y directions.

if mod(N(1),2) == 0
    qx = circshift(((-N(1)/2):(N(1)/2-1))/(N(1)*pixelSize),[0 -N(1)/2]);
else
    qx = circshift((((1-N(1))/2):((N(1)-1)/2))/(N(1)*pixelSize),[0 (1-N(1))/2]);
end
if nargout == 2
    % Check to see if second dimension length is same as first
    if length(N) == 1
        [qy,qx] = meshgrid(qx);
    else
        % Add second dimension
        if mod(N(2),2) == 0
            qy = circshift(((-N(2)/2):(N(2)/2-1))/(N(2)*pixelSize),[0 -N(2)/2]);
        else
            qy = circshift((((1-N(2))/2):((N(2)-1)/2))/(N(2)*pixelSize),[0 (1-N(2))/2]);
        end
        [qy,qx] = meshgrid(qy,qx);
    end
end
end