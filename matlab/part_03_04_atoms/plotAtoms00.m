function [] = plotAtoms00(atomData,atomLattice)

% Colin Ophus - clophus@lbl.gov - 2020 February
% APS tutorial example 02 - plotting atoms from the FePt sample

% In this example, we will plot atoms in 3D, with separation between atomic
% planes, coloring, shading (depth cueing)
% Version 01 - 3D scatter plot with some plotting tricks

% Inputs:
% atomData - [N x 4] array, where first three columns are (x,y,z) and last
%            column is the atom ID (1, 2, ...).
% atomLattice - [N x 4] array representing the basis for a fitted lattice:
%               [1 a b c] for the linear lattice equation:
%               r = r0 + a*u + b*v + c*w where (u,v,w) are 3 lattice vectors

% Other variables:
flagCenterAtoms = 0;
flagMoveAtomsToLattice = false;

% Appearance of atoms
atomSize = 40;
atomColors = [ ...
    1 0 0;
    0 0.7 1];
atomColorEdges = [0 0 0];
flagPlotAtomTint = true;
atomLinewidth = 1;

% tinting
atomColorTint = [1 1 1];
atomSizeTint = 5;
atomShiftTint = [0 -1 2];


% Shadowing
flagApplyShadows = 0;
% radiusShadow = 4;
% shadingShadow = 0.5;  % reduction of intensity for being shadowed
%vecShadow = [1 -1 1];  % direction vector for casting shadow
radiusShadow = [3 4.5];
shadingShadow = [0.7 0.85];  % reduction of intensity for being shadowed
vecShadow = [0.5 1 0.5];  % direction vector for casting shadow
% % moving tint reflection to match shadow
% atomSizeTint = 3;
% atomShiftTint = [0 2 1];

% Shading / depth cueing
flagApplyDepthShading = 0;
atomsDepthShade = [1 0.25];  % [closest_to_camera furthest_from_camera]
atomsDepthPower = 1;


% positioning and splitting of atomic planes
splitPlanesBunch = 2;   % How many planes to group together in each slice
splitPlanesShift = -1 - 0.5*0;  % additive shift of atomic planes
splitPlanesDistance = 70*0;  % multiplicative shift of atomic planes (spacing)
splitPlanesIndex = 2;  % 1 for x, 2 for y, 3 for z


% Camera stuff
cameraTarget = [0 0 0];
cameraPosition = [4 2 1]*1000;
cameraAngle = 5;

% init
numberAtoms = size(atomData,1);


% position of all atoms
atomPos = zeros(numberAtoms,3);
atomPos(:) = atomData(:,1:3);
if flagCenterAtoms == true
    for a0 = 1:3
        atomPos(:,a0) = atomPos(:,a0) - mean(atomPos(:,a0));
    end
end

% Shift atoms to best fit lattice coordinates
if flagMoveAtomsToLattice == true
    subFit = ~isinf(sum(atomLattice,2)) ...
        & ~isnan(sum(atomLattice,2));
    
    lat = atomLattice(subFit,:) \ atomPos(subFit,:);
    atomPos(subFit,:) = atomLattice(subFit,:) * lat;
end


% Measure distance of all atoms from camera
distCamera = sqrt( ...
    (atomPos(:,1) - cameraPosition(1)).^2 + ...
    (atomPos(:,2) - cameraPosition(2)).^2 + ...
    (atomPos(:,3) - cameraPosition(3)).^2);

% Split apart atomic planes to reveal interior structure
if splitPlanesDistance > 0
    atomIndex = atomLattice(:,splitPlanesIndex+1);
    atomIndex(:) = atomIndex - median(atomIndex) + splitPlanesShift;
    atomIndex(:) = round(atomIndex / splitPlanesBunch);
    
    atomPos(:,splitPlanesIndex) = atomPos(:,splitPlanesIndex) ...
        + atomIndex * splitPlanesDistance;
end


% Generate coloring for all atoms
atomTypes = unique(atomData(:,4));
atomRGB = zeros(numberAtoms,3);
for a0 = 1:length(atomTypes)
    sub = atomData(:,4) == atomTypes(a0);
    atomRGB(sub,1) = atomColors(a0,1);
    atomRGB(sub,2) = atomColors(a0,2);
    atomRGB(sub,3) = atomColors(a0,3);
end

% Apply depth shading / depth cueing
if flagApplyDepthShading == true
    scaleColor = distCamera - min(distCamera);
    scaleColor(:) = scaleColor / max(scaleColor);  % scaled from 0 to 1
    scaleColor(:) = 3*scaleColor(:).^2 - 2*scaleColor(:).^3;
    scaleColor(:) = scaleColor.^atomsDepthPower;
    scaleColor = atomsDepthShade(1) + ...
        scaleColor*(atomsDepthShade(2) - atomsDepthShade(1));  % scaled to atomsDepthShade
    atomRGB(:) = atomRGB .* scaleColor;
end


% Apply shadowing
if flagApplyShadows == true
    vecShadow = vecShadow / norm(vecShadow);
    distProj2 = zeros(numberAtoms,1);
    sub = false(numberAtoms,1);
    vecRep = repmat(vecShadow,[numberAtoms 1]);
    radiusShadow2 = radiusShadow.^2;
    
    subShadowed = false(numberAtoms,1);
    
    % Loop through each atomic site, determine if it is shadowed or not
    for a0 = 1:numberAtoms
        xyz0 = atomPos(a0,1:3);  % site under consideration
        
        % Projected distance from shadow vector direction starting from xyz0
        %         distProj(:) = abs(sum((atomPos - xyz0).*vecShadow,2));
        distProj2(:) = sum(cross(atomPos - xyz0, vecRep).^2,2);
        % subset of atoms on the correct side to be shadowed
        sub(:) = sum((atomPos - xyz0).*vecShadow,2) > 0;
        minDist2 = min(distProj2(sub));
        
        if minDist2 < radiusShadow2(1)
            atomRGB(a0,:) = atomRGB(a0,:) * shadingShadow(1);
            subShadowed(a0) = true;
        elseif minDist2 < radiusShadow2(2)
            atomRGB(a0,:) = atomRGB(a0,:) * shadingShadow(2);
            subShadowed(a0) = true;
        end
    end
end



% plotting
figure(12)
clf
% set(gcf,'outerposition',[500 500 576+512 512])
% set(gcf,'outerposition',[500 500 576+512-6 512+6])
% set(gcf,'color','w');
set(gcf,'outerposition',[500 500 576+512-5 512+8])

% axes('position',[0 0 1 1])
hold on

% Plot atoms
scatter3( ...
    atomPos(:,1),...
    atomPos(:,2),...
    atomPos(:,3));
% scatter3( ...
%     atomPos(:,1),...
%     atomPos(:,2),...
%     atomPos(:,3),...
%     ones(numberAtoms,1)*atomSize,...
%     atomRGB);
% scatter3( ...
%     atomPos(:,1),...
%     atomPos(:,2),...
%     atomPos(:,3),...
%     ones(numberAtoms,1)*atomSize,...
%     atomRGB,...
%     'filled','marker','o','linewidth',atomLinewidth,...
%     'markeredgecolor',atomColorEdges);


if flagPlotAtomTint == true
    % Move tints in slightly front of other atoms, apply shift
    shiftTint = cameraPosition / norm(cameraPosition);
    shiftTint = shiftTint*0.1 + atomShiftTint;
    
    %     if flagApplyShadows == false
    %         scatter3( ...
    %             atomPos(:,1) + shiftTint(1),...
    %             atomPos(:,2) + shiftTint(2),...
    %             atomPos(:,3) + shiftTint(3),...
    %             'sizedata',atomSizeTint,...
    %             'marker','o',...
    %             'markeredgecolor','none',...
    %             'markerfacecolor',atomColorTint);
    %     else
    %
    %        scatter3( ...
    %             atomPos(~subShadowed,1) + shiftTint(1),...
    %             atomPos(~subShadowed,2) + shiftTint(2),...
    %             atomPos(~subShadowed,3) + shiftTint(3),...
    %             'sizedata',atomSizeTint,...
    %             'marker','o',...
    %             'markeredgecolor','none',...
    %             'markerfacecolor',atomColorTint);
    %     end
end



hold off

% axis equal off


% % Camera stuff
% camtarget(cameraTarget)
% campos(cameraTarget + cameraPosition);
% camva(cameraAngle)
% camproj('perspective')




end