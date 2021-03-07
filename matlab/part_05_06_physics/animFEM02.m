function [] = animFEM02(fbase)

% Colin Ophus - clophus@lbl.gov - 2020 February
% APS tutorial example - Animation to show interaction of STEM probe with a
%                        sample that ranges from fully ordered on the left
%                        side, to fully disordered on the right side.
% 02 - This script calls 01, and animates the ordered --> disordered 
% transformation of the sample.


% Ordering animation
Nc = 60;
probeSites = [5 45];
rScaleArray = linspace(0,1,Nc) * 0.0015;

for a0 = 1:length(rScaleArray)
    rScale = rScaleArray(a0);
    
    animFEM01(probeSites,rScale);
    
    if nargin > 0
        fname = [fbase sprintf('%04d',a0) '.png'];
        eval(['export_fig -nocrop -r110 ' fname]);
    else
        drawnow;
    end
end

end