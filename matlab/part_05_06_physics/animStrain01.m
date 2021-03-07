function [] = animStrain01(probeSites)

% Colin Ophus - clophus@lbl.gov - 2020 February
% APS tutorial example - Animation to show interaction of STEM probe with a
%                        crystalline sample, where the sample has been
%                        compressed on one side, and stretched on the other.
% 01 - This script renders the scene.

skip = 1;
flagScatter = 1;


imageSize = [1 1]*128;
probeRadius = 0.6;
% probeSites = [ ...
%     5 5;
%     ...%5 15;
%     5 25;
%     ...%5 35;
%     5 45];
probeSites(:,2) = probeSites(:,2) + 0.2;
probePos = [20 0 -20]*1.2 + 5;
% probeAngle = 2;
probeXY = [-5 5;-5 5];
probeRep = 3;
probeCrop = imageSize;%[1 1]*100;
probeImageScale = 16;%8*2;
probeIntRange = [3 30];
linewidthProbes = 1;

numAtoms = [10 50 10];
cellDim = [11 51 11];

mSize = 40 + 10;
lw = 1;
radius = 0.48;

cTar = [0 0 0] + cellDim * 0.5 + [0 -1.3 0];
% cPos = [-3 0 2]*32;
cPos = [-3 0 1.5]*36;
cAngle = 20 + 1;


% Atomic sites
[yp,xp,zp] = meshgrid(1:numAtoms(2),1:numAtoms(1),1:numAtoms(3));
p = [xp(:) yp(:) zp(:)];
% Shifts
dy = (p(:,2) - 0).^2*0.008;
% dy = (p(:,2) - 0).^2*0.009 - p(:,2)*0.0;
p(:,2) = p(:,2)*0.6 + dy + 0.5;
p(:,2) = numAtoms(2)+1-p(:,2);

% Projected sites for FFTs
% x = 1:probeCrop(1);
% y = 1:probeCrop(2);
x = linspace(-probeImageScale/2,probeImageScale/2,probeCrop(1));
y = linspace(-probeImageScale/2,probeImageScale/2,probeCrop(2));
[yi,xi] = meshgrid(y-mean(y),x-mean(x));
w2 = hanning(imageSize(1))*hanning(imageSize(2))';
% w2 = sqrt(w2);
k = fspecial('gaussian',21,1.5);
kDisk = fspecial('disk',8-1);

% Sphere
% [xs,ys,zs] = sphere(8*2);
% xs = xs * radius;
% ys = ys * radius;
% zs = zs * radius;


% Plotting
h = figure(1);
clf
% set(h,'color','w','outerposition',[1020 380 576 512+128])
set(h,'color','w')
hold on
% sub = p(:,1) < 1.5 - 100;
% scatter3(p(sub,2),p(sub,1),p(sub,3),...
%     'marker','o','sizedata',mSize,'linewidth',lw,...
%     'markerfacecolor',[1 0.7 0.7],'markeredgecolor',[1 0 0])
%
% scatter3(p(~sub,2),p(~sub,1),p(~sub,3),...
%     'marker','o','sizedata',mSize,'linewidth',lw,...
%     'markerfacecolor',[1 0.9 0.9],'markeredgecolor',[1 0.5 0.5])

if flagScatter == true
%     scatter3(p(:,2),p(:,1),p(:,3),...
%         'marker','o','sizedata',mSize,'linewidth',lw,...
%         'markerfacecolor',[1 0.7 0.7],'markeredgecolor',[1 0 0])
    
    
    
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
        if sub(a0) == true
            %         c1 = [1 0.7 0.7];
            c2 = [1 0 0];
        else
            %         c1 = [1 0.9 0.9];
            c2 = [0.5 0 0];
        end
        
        surf(xs+xyz(2),ys+xyz(1),zs+xyz(3),...
            'edgecolor','none','facecolor',c2)
        
    end
end

% Draw probes
t = linspace(0,2*pi,180+1);
ct = cos(t);
st = sin(t);

% imageProbe = zeros(imageSize);
for a0 = 1:size(probeSites,1)
    xR = probeXY(1,:) + probeSites(a0,1);
    yR = probeXY(2,:) + probeSites(a0,2);
    sub = p(:,1) >= xR(1) ...
        & p(:,1) <= xR(2) ...
        & p(:,2) >= yR(1) ...
        & p(:,2) <= yR(2);
    
    xInd = (p(sub,1) - xR(1)) / (xR(2) - xR(1));
    yInd = (p(sub,2) - yR(1)) / (yR(2) - yR(1));
    xInd = min(max(round(xInd * imageSize(1)),1),imageSize(1));
    yInd = min(max(round(yInd * imageSize(2)),1),imageSize(2));
    imageProbe = accumarray([xInd yInd],...
        ones(sum(sub),1),imageSize);
    imageProbe(:) = conv2(imageProbe,k,'same');
    imageProbe(:) = imageProbe .* w2;
    imageProbe = repmat(imageProbe,[1 1]*probeRep);
    imageProbe(:) = fft2(imageProbe);
    imageProbe = fftshift(abs(imageProbe( ...
        [(1:probeCrop(1)) ((1-probeCrop(1)):0)+size(imageProbe,1)],...
        [(1:probeCrop(2)) ((1-probeCrop(2)):0)+size(imageProbe,2)])));
    imageProbe(:) = conv2(imageProbe,kDisk,'same');
    
    
    warp(yi+probeSites(a0,2),...
        xi+probeSites(a0,1),...
        xi*0+probePos(3),...
        imageProbe);
    line([-1 1]*probeRadius+probeSites(a0,2),...
        [0 0]+probeSites(a0,1),...
        probePos([1 3]),...
        'linewidth',linewidthProbes,'color',[0 0.8 1])
    line([-1 1]*probeRadius+probeSites(a0,2),...
        [0 0]+probeSites(a0,1),...
        probePos([3 1]),...
        'linewidth',linewidthProbes,'color',[0 0.8 1])
    line(ct*probeRadius+probeSites(a0,2),...
        st*probeRadius+probeSites(a0,1),...
        probePos(3)+ct*0,...
        'linewidth',linewidthProbes,'color',[0 0.8 1])
    
end
colormap(gray(256))
caxis(probeIntRange)

hold off
axis equal off

camtarget(cTar([2 1 3]))
campos(cPos([2 1 3])+cTar([2 1 3]))
camva(cAngle)
camproj('perspective')

% camlight right
% light('style','infinite','position',[0 -6 1]*1e3)
% light('style','infinite','position',[0 -6 -1]*1e3)
% lighting phong
% material shiny

% camlight left
% light('style','infinite','position',[1-15 -6 1+1 + 5]*1e3)
% light('style','infinite','position',[0-15 -6 -1+1 + 5]*1e3)


% figure(678)
% clf
% imagesc(imageProbe)
% axis equal off
% colormap(gray(256))

% toc
end