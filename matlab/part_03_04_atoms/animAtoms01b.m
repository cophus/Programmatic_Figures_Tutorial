function [] = animAtoms01b(atomData,atomLattice,fileNameBase)

% Colin Ophus - clophus@lbl.gov - 2020 February
% APS tutorial example - animating the atom plot from the FePt sample

% This script calls animAtoms01a.m with rotating camera angles, and then
% outputs images with export_fig.m in order to render a movie.

% Inputs
% fbase - Character string for output file names.  Note that you must
%         create output directories before writing to them.
numberFrames = 60*6;
thetaStart = 26.565*pi/180;
cameraRadius = norm([4 2])*1000;
cameraHeight = 1000;


% camera variables
thetaArray = linspace(0,2*pi,numberFrames+1) + thetaStart;
thetaArray(end) = [];

% Main loop
for a0 = 1:numberFrames
    cameraPosition = [ ...
        cameraRadius*cos(thetaArray(a0)) ...
        cameraRadius*sin(thetaArray(a0)) ...
        cameraHeight];
    
    % draw figure
    animAtoms01a(atomData,atomLattice,cameraPosition);
    
    % Output image
    fileName = [fileNameBase sprintf('%04d',a0) '.png'];
    evalc(['export_fig -nocrop -r90 ' fileName]);
end



end