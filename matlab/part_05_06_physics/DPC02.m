function [sDPC] = DPC02(sDPC,fbase)

% Colin Ophus - clophus@lbl.gov - 2020 February
% APS tutorial example - DPC plot / movies

flagPlot = 1;

% convergence semiangle - from 0.02 to 0.32 roughly
sDPC.qProbeMax = 0.16;


sDPC.probeRange = sDPC.imageSize(1) * [0 1];
sDPC.numProbes = 60*4*4 / 4;
% vecInds = 480+8;
vecInds = 1:sDPC.numProbes;

% sDPC.qProbeMax = 0.32;
% sDPC.qProbeMax = 0.16;
% sDPC.qProbeMax = 0.06;
% sDPC.qProbeMax = 0.02;

yR1 = [0 pi]*1.1;
xCut = 512;
yCut = 182;

sDPC.sigmaPot = 1;
sigCountsTotal = 500*0;%1e3;


intRangeFFT = [0.0 1.1];
intPowerFFT = 1;
intRange = [0 0.004];
intPower = 1;


% scaling
intRange = intRange / (0.06 / sDPC.qProbeMax)^2;

% Probes
sDPC.xp = linspace(sDPC.probeRange(1),sDPC.probeRange(2),sDPC.numProbes+1);
sDPC.xp(end) = [];

% Potentials
trans = exp((1i*sDPC.sigmaPot)*sDPC.pot);


% Generate initial probe
sDPC.q2 = sDPC.qxa.^2 + sDPC.qya.^2;
sDPC.q1 = sqrt(sDPC.q2);
dq = sDPC.qxa(2,1);
sDPC.psi0 = min(max((sDPC.qProbeMax - sDPC.q1)/dq+0.5,0),1);

% signal counts
sigCounts = sigCountsTotal / sum(abs(sDPC.psi0(:).^2));


% Coords
xc = [(1:(xCut/2)) ((1-xCut/2):0)+sDPC.imageSize(1)];
yc = [(1:(yCut/2)) ((1-yCut/2):0)+sDPC.imageSize(2)];



% Main
% Psi0 = fft2(sDPC.psi0);
psi0 = zeros(sDPC.imageSize);
psi = zeros(sDPC.imageSize);
Psi = zeros(sDPC.imageSize);

% Main loop

sDPC.CoM = zeros(sDPC.numProbes,2);

for a0 = vecInds
    xx = sDPC.xp(a0);
    
    
    qShift = exp((-2i*pi*xx)*sDPC.qxa);
    psi0(:) = ifft2(sDPC.psi0 .* qShift);
    psi(:) = psi0 .* trans;
    Psi(:) = fft2(psi);
    
    PsiInt = abs(Psi).^2;
    if sigCounts > 0
        PsiInt(:) = poissrnd(PsiInt*sigCounts)/sigCounts;
    end
    xCoM = sum(PsiInt(:).*sDPC.qxa(:)) / sum(PsiInt(:));
    sDPC.CoM(a0,:) = [sDPC.xp(a0) xCoM];
    
    if flagPlot == true
        P = exp(1i*pi) * Psi ./ qShift;
        Pc = fftshift(P(xc,yc));
        
        p = psi(:,1);
        pamp = abs(p).*intPower;
        prgb = colorComplex(ones(size(pamp)),angle(p),[0 1]);
        pamp(:) = (pamp - intRange(1).^intPower) ...
            / (intRange(2).^intPower - intRange(1).^intPower);
        pamp(:) = min(max(pamp,0),pi);
        
        if sigCounts > 0
            Pint = abs(Pc).^2;
            Pint(:) = poissrnd(Pint*sigCounts)/sigCounts;
            Pc(:) = sqrt(Pint).*exp(1i*angle(Pc));
        end
        
        Igray = abs(Pc).^intPowerFFT;
        Igray(:) = (Igray - intRangeFFT(1).^intPowerFFT) ...
            / (intRangeFFT(2).^intPowerFFT - intRangeFFT(1).^intPowerFFT);
        Igray = repmat(Igray,[1 1 3]); 
        Irgb = colorComplex(abs(Pc).^intPowerFFT,...
            angle(Pc),intRangeFFT.^intPowerFFT);
        
        
        figure(11)
        clf
        set(gcf,'color','w')
        %         set(gcf,...
        %             'outerposition',[1000 500-256 576 512+256],...
        %             'color','w')
        
        b = [0.01 0.004];
        axes('position',[0+b(1) 0.6+b(2) 1-2*b(1) 0.4-2*b(2)])
        %         plot(sDPC.x,sDPC.p,'linewidth',2,'color','k')
        area(sDPC.x(:),...
            sDPC.p(:),...
            'linewidth',2,...
            'facecolor',[1 1 1]*0.6,'edgecolor','none')
        hold on
        %         area(sDPC.x,pamp,'linewidth',2,...
        %             'facecolor',prgb,'edgecolor','none')
        % %         patch([sDPC.x(:); sDPC.x(end); sDPC.x(1)],...
        % %             [pamp(:); 0; 0],'r',...
        % %             'facecolor','interp',...
        % %             'facevertexcdata',...
        % %             [squeeze(prgb); 0 0 0; 0 0 0])
        %         xp = [sDPC.x(:); flipud(sDPC.x(:))];
        %         yp = [pamp(:); pamp(:)*0];
        v = [sDPC.x(:) pamp(:);sDPC.x(:) pamp(:)*0];
        Nv = size(v,1)/2;
        f = [(1:Nv-1)' (1:Nv-1)'+1 ...
            (1:Nv-1)'+1+Nv (1:Nv-1)'+Nv];
        patch('vertices',v,'faces',f,...
            'facecolor','interp',...
            'edgecolor','none',...
            'cdata',repmat(prgb,[1 2 1]))
        
        hold off
        
        box on
        set(gca,'xtick',[],'ytick',[])
        xlim(sDPC.x([1 end]))
        ylim(yR1)
        
        
        axes('position',[0+b(1) 0.3+b(2) 1-2*b(1) 0.3-2*b(2)])
        imagesc(permute(Igray,[2 1 3]))
        hold on
        xx = xCoM*sDPC.imageSize(1);
        plot([0 0]+xx+sDPC.imageSize(1)/2+1,...
            [1 yCut],...
            'linewidth',1,'color',[0 0.8 0])
        
        hold off
        axis equal off
        
        
        axes('position',[0+b(1) 0.0+b(2) 1-2*b(1) 0.3-2*b(2)])
        imagesc(permute(Irgb,[2 1 3]))
        axis equal off
        
        drawnow
        
        if nargin > 1
            fname = [fbase sprintf('%04d',a0) '.png'];
            eval(['export_fig -nocrop -r101 ' fname]);
        end
        
%         %         Ip = abs(Psi).^intPowerFFT;
%         P = fftshift(exp(1i*pi) * Psi ./ qShift);
%         Irgb = colorComplex(abs(P).^intPowerFFT,...
%             angle(P),intRangeFFT.^intPowerFFT);
%         imagesc(Irgb)
%         %imagesc(fftshift(Ip.^intPowerFFT));
%         %         imagesc([abs(psi0) abs(psi)])
%         %         imagesc(fftshift(sDPC.psi0))
%         axis equal off
%         drawnow
%         %colormap(gray(256))
%         %         caxis(intRangeFFT.^intPowerFFT)
        

    else
        if length(vecInds) > 1
            comp = a0 / sDPC.numProbes;
            progressbar(comp,2);
        end
    end
    
        
    
end


end