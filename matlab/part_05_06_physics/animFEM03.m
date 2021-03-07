function [] = animFEM03(fbase)

% Colin Ophus - clophus@lbl.gov - 2020 February
% APS tutorial example - Animation to show interaction of STEM probe with a
%                        sample that ranges from fully ordered on the left
%                        side, to fully disordered on the right side.
% 03 - This script calls 01, animates STEM probe motion.


% Ordering animation
Nc = 120;
v = linspace(45,5,Nc+1)';
probeSitesArray = [ ...
    [ones(Nc,1)*5 v(1:(end-1))];
    [ones(Nc,1)*5 flipud(v(2:end))]
    ];
rScale = 0.0015;

for a0 = 1:size(probeSitesArray,1)
    probeSites = probeSitesArray(a0,:);
    animFEM01(probeSites,rScale);
    
    if nargin > 0
        fname = [fbase sprintf('%04d',a0) '.png'];
        eval(['export_fig -nocrop -r110 ' fname]);
    else
        drawnow;
    end
end

end