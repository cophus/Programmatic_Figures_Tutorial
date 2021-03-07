function [sDPC] = DPC01()

% Colin Ophus - clophus@lbl.gov - 2020 February
% APS tutorial example - DPC plot / movies

% simulate a CoM DPC reconstruction as a 1D slice

sDPC.pixelSize = 1;
sDPC.imageSize = [1 1]*512*1;

% Coordinates
[sDPC.qxa,sDPC.qya] = makeFourierCoords(sDPC.imageSize,sDPC.pixelSize);



% Generate atomic potential
x = (1:sDPC.imageSize(1))';

p = zeros(sDPC.imageSize(1),1);

numPeaks = 20;
t = linspace(0,1,numPeaks);
% xp = sDPC.imageSize(1)*0.34 + t*sDPC.imageSize(1)*0.6;
% Ip = t*pi;
% sigma = sDPC.imageSize(1)*0.009;
xp = sDPC.imageSize(1)*0.25 + t*sDPC.imageSize(1)*0.5;
Ip = sqrt(1-(1-2*t).^2)*pi;
sigma = sDPC.imageSize(1)*0.008;

% sigma = sDPC.imageSize(1)*0.002;


for a0 = 1:length(t)
    p(:) = p + ...
        Ip(a0)*exp(-(x-xp(a0)).^2/(2*sigma^2));
end
sDPC.x = x;
sDPC.p = p;
sDPC.pot = repmat(p,[1 sDPC.imageSize(2)]);



figure(11)
clf
plot(x,p,'linewidth',2,'color','r')



end