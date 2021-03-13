function [] = animOrient01(probeSites)

% Colin Ophus - clophus@lbl.gov - 2020 February
% APS tutorial example - Animation to show interaction of STEM probe with a
%                        crystalline sample, consisting of three grains 
%                        with different orientations.
% 01 - This script renders the scene.

rng(1);  % set random see for repeatable behaviour

skip = 1;
flagScatter = 1;

imageSizeOutput = [200 200];
numAtoms = [10 50 10];
pixelSize = 0.05;
sigmaAtoms = 4;

linewidthProbes = 1.5;%1.2 - 0.3;
probeRadius = 0.7;
if nargin == 0
probeSites = [ ...
    5 05;
    ... % 5 15;
    5 25;
    ... % 5 35;
    5 45];
end
% probeWeights = [1.0 0.75 0.50 0.25 0];
probeRep = 3;

% appearance
mSize = 50;% * 4;
mSizeTint = 10;%*1.5;
lw = 1;
radius = 0.48 + 0.04;
% radius = 0.36;

probePos = [20 0 -20]*1.2 + 5;
probeImageScale = 15;

% Camera
cellDim = numAtoms + 0;
cTar = [0 0 0] + cellDim * 0.5 + [0 -0.6 0];
cPos = [-3 0 1.5]*36;
cAngle = 20 + 1;
% % More dramatic camera angles
% cellDim = numAtoms + 0;
% cTar = [0 0 0] + cellDim * 0.5 + [10 5 -15];
% cPos = [-2.25 1.75 0]*13;
% cAngle = 75;


% Atomic sites
numAtomsMax = max(numAtoms);
v = -round(numAtomsMax/2):round(numAtomsMax/2);
[yy,xx,zz] = meshgrid(v,v,v);
p0 = [xx(:) yy(:) zz(:)];
Np = size(p0,1);

p1 = p0;
p2 = p0;
p3 = p0;

t = [45 54 -18.43]*pi/180;
m1 = [cos(t(1)) -sin(t(1)) 0;
    sin(t(1)) cos(t(1)) 0;
    0 0 1];
m2 = [1 0 0;
    0 cos(t(2)) -sin(t(2));
    0 sin(t(2)) cos(t(2))];
m3 = [cos(t(3)) -sin(t(3)) 0;
    sin(t(3)) cos(t(3)) 0;
    0 0 1];
p2 = p2 * m1 * m2 * m3;


% t = [32 47 -70]*pi/180;
t = [64 45+4-1 -70]*pi/180;
m1 = [cos(t(1)) -sin(t(1)) 0;
    sin(t(1)) cos(t(1)) 0;
    0 0 1];
m2 = [1 0 0;
    0 cos(t(2)) -sin(t(2));
    0 sin(t(2)) cos(t(2))];
m3 = [cos(t(3)) -sin(t(3)) 0;
    sin(t(3)) cos(t(3)) 0;
    0 0 1];
p1 = p1 * m1 * m2 * m3;

    

p1 = p1 + repmat([0 0 0],[Np 1]);
p2 = p2 + repmat([0 20 0],[Np 1]);
p3 = p3 + repmat([0 40 0],[Np 1]);

% abcd12 = [0.2 1 1.2 22];
% abcd23 = [0.8 1 -0.6 37];


abcd12 = [0.2 1 1.2 / 2 22];
abcd23 = [0.8 1 -0.6 / 2 37];

del = p1(:,1) < 1 | p1(:,1) > cellDim(1) ...
    | p1(:,2) < 1 | p1(:,2) > cellDim(2) ...
    | p1(:,3) < 1 | p1(:,3) > cellDim(3) ...
    | p1(:,1)*abcd12(1) + p1(:,2)*abcd12(2) + p1(:,3)*abcd12(3) > abcd12(4) - 0.1;
p1(del,:) = [];

del = p2(:,1) < 1 | p2(:,1) > cellDim(1) ...
    | p2(:,2) < 1 | p2(:,2) > cellDim(2) ...
    | p2(:,3) < 1 | p2(:,3) > cellDim(3) ...
    | p2(:,1)*abcd12(1) + p2(:,2)*abcd12(2) + p2(:,3)*abcd12(3) < abcd12(4) + 0.1 ...
    | p2(:,1)*abcd23(1) + p2(:,2)*abcd23(2) + p2(:,3)*abcd23(3) > abcd23(4) - 0.1;
p2(del,:) = [];

del = p3(:,1) < 1 | p3(:,1) > cellDim(1) ...
    | p3(:,2) < 1 | p3(:,2) > cellDim(2) ...
    | p3(:,3) < 1 | p3(:,3) > cellDim(3) ...
    | p3(:,1)*abcd23(1) + p3(:,2)*abcd23(2) + p3(:,3)*abcd23(3) < abcd23(4) + 0.1;
p3(del,:) = [];


p = [p1;p2; p3];
% p = p1;
% [yp,xp,zp] = meshgrid(1:numAtoms(2),1:numAtoms(1),1:numAtoms(3));
% p = [xp(:) yp(:) zp(:)];
% p = [1 1 1];


% Potential image
k = fspecial('gaussian',2*ceil(4*sigmaAtoms)+1,sigmaAtoms);
imageSize = round(cellDim(1:2) ./ pixelSize);
xInd = (p(:,1) - 0.5) / pixelSize;
yInd = (p(:,2) - 0.5) / pixelSize;
xInd = min(max(round(xInd),1),imageSize(1));
yInd = min(max(round(yInd),1),imageSize(2));
imagePot = accumarray([xInd yInd],...
    ones(length(xInd),1),imageSize);
imagePot(:) = conv2(imagePot,k,'same');

% Image coordinates
x = linspace(-probeImageScale/2,probeImageScale/2,imageSizeOutput(1));
y = linspace(-probeImageScale/2,probeImageScale/2,imageSizeOutput(2));
[yi,xi] = meshgrid(y-mean(y),x-mean(x));


% Make probes
[qxa,qya] = makeFourierCoords(imageSize,1);
q2 = qxa.^2 + qya.^2;
qMask = q2 < 0.25^2;


% Other image stuff
% w2 = hanning(imageSizeOutput(1)) * hanning(imageSizeOutput(2))';
w2 = tukeywinMake(imageSizeOutput);
% xx = 1:imageSizeOutput(1);
% yy = 1:imageSizeOutput(2);
% [ya,xa] = meshgrid(
xx = linspace(-1,1,imageSizeOutput(1));
yy = linspace(-1,1,imageSizeOutput(2));
[ya,xa] = meshgrid(yy,xx);
ra2 = xa.^2 + ya.^2;
ra = sqrt(ra2);
% wa = max(1 - ra,0);
% wa = exp(-ra2);
kDisk = fspecial('disk',8-2);
wa = 1 ...
    + 10*ra.*exp(-(ra - 0.35).^2/(2*0.07^2)) ...
    + 12*ra.*exp(-(ra - 0.70).^2/(2*0.09^2)) ...
    + 40*ra.*exp(-(ra - 1.05).^2/(2*0.13^2));

% Sphere
[xs,ys,zs] = sphere(8);
xs = xs * radius;
ys = ys * radius;
zs = zs * radius;



% Plotting
h = figure(1);
clf
set(h,'color','w')
% set(h,'color','w','outerposition',[1020 380-120 576 512+128])
% set(h,'color','w','outerposition',[100 100 576+576 512+512])
hold on
% sub = p(:,1) < 1.5 - 100*1;
% Np = size(p,1);
if flagScatter == true
    scatter3(p3(:,2),p3(:,1),p3(:,3),...
        'marker','o','sizedata',mSize,'linewidth',lw,...
        'markeredgecolor',[1 0.7 0.7]*0,'markerfacecolor',[1 0 0])
    scatter3(p2(:,2),p2(:,1),p2(:,3),...
        'marker','o','sizedata',mSize,'linewidth',lw,...
        'markeredgecolor',[0.8 1 0.8]*0,'markerfacecolor',[0 0.8 0])
    scatter3(p1(:,2),p1(:,1),p1(:,3),...
        'marker','o','sizedata',mSize,'linewidth',lw,...
        'markeredgecolor',[1.0 0.5 1]*0,'markerfacecolor',[1 0 1])
    
    
    d = cPos;
    d = d / norm(d) * 0.05 + [0 -0.1 0.1];
    
    scatter3(p1(:,2)+d(1),p1(:,1)+d(2),p1(:,3)+d(3),...
        'marker','o','sizedata',mSizeTint,'linewidth',lw,...
        'markerfacecolor',[1 0.8 1],'markeredgecolor','none')
    scatter3(p2(:,2)+d(1),p2(:,1)+d(2),p2(:,3)+d(3),...
        'marker','o','sizedata',mSizeTint,'linewidth',lw,...
        'markerfacecolor',[0.8 1 0.8],'markeredgecolor','none')
    scatter3(p3(:,2)+d(1),p3(:,1)+d(2),p3(:,3)+d(3),...
        'marker','o','sizedata',mSizeTint,'linewidth',lw,...
        'markerfacecolor',[1 0.9 0.9],'markeredgecolor','none')
    
    
else
    c = [0 0.4 0.6];
    for a0 = 1:skip:size(p1,1)
        xyz = p1(a0,:);
        surf(xs+xyz(2),ys+xyz(1),zs+xyz(3),...
            'edgecolor','none','facecolor',c)
    end
    
    c = [0 0.3 0];
    for a0 = 1:skip:size(p2,1)
        xyz = p2(a0,:);
        surf(xs+xyz(2),ys+xyz(1),zs+xyz(3),...
            'edgecolor','none','facecolor',c)
    end
    
    c = [0.3 0 0];
    for a0 = 1:skip:size(p3,1)
        xyz = p3(a0,:);
        surf(xs+xyz(2),ys+xyz(1),zs+xyz(3),...
            'edgecolor','none','facecolor',c)
    end
end


% Probes coords
t = linspace(0,2*pi,180+1);
ct = cos(t);
st = sin(t);
Nprobes = size(probeSites,1);
xVec = [(1:(imageSizeOutput(1)/2)) ...
    ((1-imageSizeOutput(1)/2):0) + imageSizeOutput(1)*probeRep];
yVec = [(1:(imageSizeOutput(2)/2)) ...
    ((1-imageSizeOutput(2)/2):0) + imageSizeOutput(2)*probeRep];


for a0 = 1:Nprobes
    line([-1 1]*probeRadius+probeSites(a0,2),...
        [0 0]+probeSites(a0,1),...
        probePos([1 3]),...
        'linewidth',linewidthProbes,'color',[0 0.8 1])
    line([1 -1]*probeRadius+probeSites(a0,2),...
        [0 0]+probeSites(a0,1),...
        probePos([1 3]),...
        'linewidth',linewidthProbes,'color',[0 0.8 1])
    line(ct*probeRadius+probeSites(a0,2),...
        st*probeRadius+probeSites(a0,1),...
        probePos(3)+ct*0,...
        'linewidth',linewidthProbes,'color',[0 0.8 1])
    
    xInds = mod((1:imageSizeOutput(1)) ...
        - imageSizeOutput(1)/2 ...
        + probeSites(a0,1)/pixelSize-1, imageSize(1))+1;
    yInds = mod((1:imageSizeOutput(2)) ...
        - imageSizeOutput(2)/2 ...
        + probeSites(a0,2)/pixelSize-1, imageSize(2))+1;
    imageCut = imagePot(round(xInds),round(yInds));
%     imageCut = fft2(imageCut .* w2);
%     imageCut = fftshift(wa .* abs(imageCut).^1);
    imageCut = fft2(repmat(imageCut .* w2,[1 1]*probeRep));
    imageCut = imageCut(xVec,yVec);
    imageCut = fftshift(abs(imageCut));
    imageCut(:) = conv2(imageCut,kDisk,'same');
    imageCut(:) = imageCut.^(1/2);
%     imageCut(:) = (imageCut).^0.667 ...
%         .* (wa * probeWeights(a0) + (1-probeWeights(a0)));
    
    warp(yi+probeSites(a0,2),...
        xi+probeSites(a0,1),...
        xi*0+probePos(3),...
        imageCut * 0.08 * 1)
end

hold off
axis equal off
caxis([0 0.25])



camlight left
light('style','infinite','position',[1 -6 1+1]*1e3)
light('style','infinite','position',[0 -6 -1+1]*1e3)
lighting phong
material shiny

camva(cAngle)
camproj('perspective')
camtarget(cTar([2 1 3]))
campos(cPos([2 1 3])+cTar([2 1 3]))


% toc
end