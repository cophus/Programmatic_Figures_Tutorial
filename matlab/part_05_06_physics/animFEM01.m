function [] = animFEM01(probeSites,rScale)

% Colin Ophus - clophus@lbl.gov - 2020 February
% APS tutorial example - Animation to show interaction of STEM probe with a
%                        sample that ranges from fully ordered on the left
%                        side, to fully disordered on the right side.
% 01 - This script renders the scene.

rng(1+2);

skip = 1;
flagScatter = 1;

padding = 50;
imageSizeOutput = [200 200];
numAtoms = [10 50 10];
pixelSize = 0.05;
sigmaAtoms = 4;

linewidthProbes = 1.5;
probeRadius = 0.3 + 0.05;
% probeSites = [ ...
%     5 05;
%     5 15;
%     5 25;
%     5 35;
%     5 45];
% probeWeights = [1.0 0.75 0.50 0.25 0];
probeRep = 3;

% appearance
mSize = 50;
lw = 1;
radius = 0.48 + 0.04;
% radius = 4;


probePos = [20 0 -20]*1.2 + 5;
probeImageScale = 10 + 5;

% Camera
cellDim = numAtoms + 0;
cTar = [0 0 0] + cellDim * 0.5 + [0 -0.6 0];
cPos = [-3 0 1.5]*36;
cAngle = 20 + 1;


% Atomic sites
[yp,xp,zp] = meshgrid(1:numAtoms(2),1:numAtoms(1),1:numAtoms(3));
p = [xp(:) yp(:) zp(:)];
for a0 = 1:size(p,1)
    if nargin == 1
        r = (numAtoms(2) - p(a0,2)).^1.5*0.0015;
    else
        r = (numAtoms(2) - p(a0,2)).^1.5*rScale;        
    end
    dxy = randn(1,2) * r;
    p(a0,1:2) = p(a0,1:2) + dxy;
end


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


% Other image stuff
w2 = hanning(imageSizeOutput(1)) * hanning(imageSizeOutput(2))';
xx = linspace(-1,1,imageSizeOutput(1));
yy = linspace(-1,1,imageSizeOutput(2));
[ya,xa] = meshgrid(yy,xx);
ra2 = xa.^2 + ya.^2;
ra = sqrt(ra2);
kDisk = fspecial('disk',2);
wa = 1 ...
    + 10*ra.*exp(-(ra - 0.35).^2/(2*0.07^2)) ...
    + 12*ra.*exp(-(ra - 0.70).^2/(2*0.09^2)) ...
    + 40*ra.*exp(-(ra - 1.05).^2/(2*0.13^2));
% wa(:) = wa * 1.2;
wa(:) = wa * 0.8;

% Sphere
[xs,ys,zs] = sphere(8);
xs = xs * radius;
ys = ys * radius;
zs = zs * radius;

% Plotting
h = figure(1);
clf
set(h,'color','w')
% set(h,'color',[1 1 1],'outerposition',[1020 380-120 576 512+128])
% set(h,'color',[1 1 1],'outerposition',[1020 260 576-2 640+10])
hold on
% sub = p(:,1) < 1.5 - 100*1;
if flagScatter == true
    %     scatter3(p(sub,2),p(sub,1),p(sub,3),...
    %         'marker','o','sizedata',mSize,'linewidth',lw,...
    %         'markerfacecolor',[1 0.7 0.7],'markeredgecolor',[1 0 0])
    %
    %     scatter3(p(~sub,2),p(~sub,1),p(~sub,3),...
    %         'marker','o','sizedata',mSize,'linewidth',lw,...
    %         'markerfacecolor',[1 0.9 0.9],'markeredgecolor',[1 0.5 0.5])
    
    %     scatter3(p(sub,2),p(sub,1),p(sub,3),...
    %         'marker','o','sizedata',mSize,'linewidth',lw,...
    %         'markerfacecolor',[1 0 0],'markeredgecolor',[0 0 0])
    
    scatter3(p(:,2),p(:,1),p(:,3),...
        'marker','o','sizedata',mSize,'linewidth',lw,...
        'markerfacecolor',[1 0 0],'markeredgecolor',[0 0 0])
    
    d = cPos;
    d = d / norm(d) * 0.05 + [0 -0.1 0.1];
    scatter3(p(:,2)+d(1),p(:,1)+d(2),p(:,3)+d(3),...
        'marker','o','sizedata',10,'linewidth',lw,...
        'markerfacecolor',[1 0.7 0.7],'markeredgecolor','none')
    
    d = cPos;
    d = d / norm(d) * 0.1 + [0 -0.1 0.1];
    scatter3(p(:,2)+d(1),p(:,1)+d(2),p(:,3)+d(3),...
        'marker','o','sizedata',3,'linewidth',lw,...
        'markerfacecolor',[1 0.9 0.9],'markeredgecolor','none')
else
    Np = size(p,1);
    for a0 = 1:skip:Np
        xyz = p(a0,:);
        c2 = [0.3 0 0];
        surf(xs+xyz(2),ys+xyz(1),zs+xyz(3),...
            'edgecolor','none','facecolor',c2)
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
    %     imageCut = imagePot(xInds,yInds);
    imageCut = imagePot(round(xInds),round(yInds));

    
    imageCut = fft2(repmat(imageCut .* w2,[1 1]*probeRep));
    imageCut = imageCut(xVec,yVec);
    imageCut = fftshift(abs(imageCut));
    imageCut(:) = conv2(imageCut,kDisk,'same');
    %     imageCut(:) = (imageCut).^0.667 ...
    %         .* (wa * probeWeights(a0) + (1-probeWeights(a0)));
    probeWeightScan = (probeSites(2) - 5) / 40;
    probeWeightScan = 1 - 0.5*min(max(probeWeightScan,0),1);
    imageCut(:) = 2*(imageCut).^0.5 ...
        .* (wa * probeWeightScan + (1-probeWeightScan));
    
    imageCut(:) = imageCut * 0.06;
    imageCut(:) = min(max(imageCut,0),1);
    %     imageCut(:) = 1 -imageCut;
    imagePad = padarray(imageCut,[1 1]*padding/2,0,'both');
    
    warp(yi+probeSites(a0,2),...
        xi+probeSites(a0,1),...
        xi*0+probePos(3),...
        imagePad)
end

hold off
axis equal off
caxis([0 1])
set(gca,'ydir','reverse')


% camlight left
% light('style','infinite','position',[1 -6 1+1]*1e3)
% light('style','infinite','position',[0 -6 -1+1]*1e3)
% lighting phong
% material shiny

camva(cAngle)
camproj('perspective')
camtarget(cTar([2 1 3]))
campos(cPos([2 1 3])+cTar([2 1 3]))

% toc
end