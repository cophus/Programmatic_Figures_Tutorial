function [EW] = makeSpiralEW(numOAM)

% Make a spiral phase EW image

imageSize = [1 1] * 1000;
qRange = [-1 1] * 0.01;
if nargin == 0
    numOAM = 0;
end


[qxa,qya] = makeFourierCoords(imageSize,1);
q2 = qxa.^2 + qya.^2;
q1 = sqrt(q2);
qt = atan2(qya,qxa);
dq = qxa(2,1) - qxa(1,1);

% Make EW
PsiAmp1 = 1 - min(max( (qRange(1) - q1)/dq + 0.5,0),1);
PsiAmp2 = min(max( (qRange(2) - q1)/dq + 0.5,0),1);
PsiAmp = PsiAmp1 .* PsiAmp2;
PsiPhase = numOAM*qt;
Psi = PsiAmp .* exp(1i*PsiPhase);
EW = fftshift(ifft2(Psi));
EW = EW / max(abs(EW(:)));

% figure(11)
% clf
% imagesc(abs(EW))
% % imagesc(angle(EW))
% % imagesc(fftshift((PsiPhase)))
% % imagesc(fftshift((PsiAmp)))
% axis equal off
% colorbar






end