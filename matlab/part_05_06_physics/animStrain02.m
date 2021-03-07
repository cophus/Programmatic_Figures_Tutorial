function [] = animStrain02(fbase)

% Colin Ophus - clophus@lbl.gov - 2020 February
% APS tutorial example - Animation to show interaction of STEM probe with a
%                        crystalline sample, where the sample has been
%                        compressed on one side, and stretched on the other.
% 02 - This script repeatedly calls 01, with different probe positions.

% Ordering animation
Nc = 120;
v = linspace(45,5,Nc+1)';
probeSitesArray = [ ...
    [ones(Nc,1)*5 v(1:(end-1))];
    [ones(Nc,1)*5 flipud(v(2:end))]
    ];
% rScale = 0.0015;

for a0 = 1:size(probeSitesArray,1)
    probeSites = probeSitesArray(a0,:);
    animStrain01(probeSites);
    
    if nargin > 0
        fname = [fbase sprintf('%04d',a0) '.png'];
        eval(['export_fig -nocrop -r110 ' fname]);
    else
        drawnow;
    end
end

end