function [] = schematic01()

% Colin Ophus - clophus@lbl.gov - 2020 February
% APS tutorial example 02 - plotting atoms from the FePt sample

zBox = [1.5 0.1 -0.1 -1.5];
% zMid = linspace(-1.5,1.5,6);
% zMid = zMid*0.7 + 0.1*zMid.^3;
% zMid = zMid*0.9 + 0.08*zMid.^3;
zMid = (linspace(0,1,5).^0.9)*3-1.5;
mxy = [0 -2.5];

cTar = [0 0 0] + [0 -1.0 0];
cPos = [-2 -0.5 1]*6;
cAngle = 15;


ampRange = [0 0.5e-4];

imageSize = [1 1]*512*2;
pixelSize = 0.02;
[qxa,qya] = makeFourierCoords(imageSize,pixelSize);
q2 = qxa.^2 + qya.^2;
prop = exp(-1i*pi*q2*0.4);
qMax = 0.8;
Psi = double(q2 <= qMax^2);
s = 0.8;
pot = (cos(qxa*s).^4) .* (cos(qya*s).^4) ...
    + (sin(qxa*s/2).^4) .* (sin(qya*s/2).^4);
pot = circshift(pot,[1 1]*16);
trans = exp(1i*pot*0.5*2);


% Testing plotting of the potential
% figure(56)
% clf
% imagesc(fftshift(pot))
% axis equal off
% colormap(gray(256))


p = [ ...
    -1 -1 0;
    1 -1 0;
    1 1 0;
    -1 1 0];
f = 1:4;



figure(1)
clf
% Uncomment these lines to fix the position and colour of the axes
% set(gcf,'outerposition',[1000 500 576 512],'color','w')
% set(gcf,'outerposition',[1000 500 576 512],'color','k')
hold on

for a0 = 1:length(zBox)
    patch('vertices',p(:,[2 1 3]).*[1 1 1] + [0 0 zBox(a0)],...
        'faces',f,'linewidth',1,...
        'facecolor',[1 1 1],'edgecolor',[0 0 0]);
end

x = linspace(0,1,imageSize(1)/2);
y = linspace(0,1,imageSize(2)/2);
[ya,xa] = meshgrid(y,x);
v = [(1:(imageSize(1)/4)) ((1-imageSize(1)/4):0)+imageSize(1)];
vec = fliplr(1:length(zMid));
scale = linspace(1,1.8,length(vec));
% bb = 8;
for a0 = 1:length(vec)
    ind = vec(a0);
    
    xyz = p(:,[2 1 3]).*[0.5 0.5 1] + [mxy zMid(ind)];
    %     patch('vertices',xyz,...
    %         'faces',f,'linewidth',1,...
    %         'facecolor',[1 1 1],'edgecolor',[0 0 0]);
    
    yy = [xyz(1,1) xyz(3,1)-xyz(1,1)];
    xx = [xyz(1,2) xyz(3,2)-xyz(1,2)];
    
    psi = ifft2(Psi);
    psi = fftshift(psi(v,v));
    Irgb = colorComplex(abs(psi),-angle(psi),ampRange);
    %     Irgb(1:bb,:,:) = 1;
    %     Irgb(:,1:bb,:) = 1;
    %     Irgb(((1-bb):0)+end,:,:) = 1;
    %     Irgb(:,((1-bb):0)+end,:) = 1;
    
    h = abs(psi).^0.8 * 100 * scale(a0);
    
    warp(ya*yy(2)+yy(1), ...
        xa*xx(2)+xx(1),...
        xa*0+xyz(1,3) + h,Irgb);
    
    Psi = Psi.*prop;
    Psi = fft2(ifft2(Psi) .* trans);
end

hold off
axis equal off
box on

camtarget(cTar)
campos(cPos + cTar)
camva(cAngle)
camproj('perspective')
set(gca,'ydir','normal')

caxis([0 1])
camlight left
material shiny




end