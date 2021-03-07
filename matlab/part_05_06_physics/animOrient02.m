function [] = animOrient02(fbase)

% Colin Ophus - clophus@lbl.gov - 2020 February
% APS tutorial example - Animation to show interaction of STEM probe with a
%                        crystalline sample, consisting of three grains 
%                        with different orientations.
% 02 - This script repeatedly calls 01, with different probe positions.

% Ordering animation
Nc = 120;   % number of probe positions / frames in each direction
v = linspace(45,5,Nc+1)';
probeSitesArray = [ ...
    [ones(Nc,1)*5 v(1:(end-1))];
    [ones(Nc,1)*5 flipud(v(2:end))]
    ];

for a0 = 1:size(probeSitesArray,1)
    probeSites = probeSitesArray(a0,:);
    animOrient01(probeSites);
    
    if nargin > 0
        fname = [fbase sprintf('%04d',a0) '.png'];
        eval(['export_fig -nocrop -r110 ' fname]);
    else
        drawnow;
    end
end

end