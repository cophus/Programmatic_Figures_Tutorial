function [] = animAtoms02b(atomData,atomLattice,fileNameBase)

% Colin Ophus - clophus@lbl.gov - 2020 February
% APS tutorial example - animating the atom plot from the FePt sample

% This script calls animAtoms01a.m with rotating camera angles, and then
% outputs images with export_fig.m in order to render a movie.

% Inputs
% fbase - Character string for output file names.  Note that you must
%         create output directories before writing to them.
numberFrames = 60*2;
splitPlanes = linspace(0,70,numberFrames);

% Main loop
for a0 = 1:numberFrames
    % draw figure
    animAtoms02a(atomData,atomLattice,splitPlanes(a0));
    
    % Output image
    fileName = [fileNameBase sprintf('%04d',a0) '.png'];
    evalc(['export_fig -nocrop -r90 ' fileName]);
end



end