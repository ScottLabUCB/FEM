function [sFEM] = FEM44(sFEM)
tic

% Colin Ophus - 2019 Sept
% 44 - elliptical polar transformation of all images

%pixelSize = sFEM.scale;
pixelSize = 0.0335; %TitanX in inverse nm
rSigma = 0.1;
tSigma = 1;
rMax = 240;
dr = 2; %1  
dt = 5; %5
sFEM.polarRadius = (0:dr:rMax)*pixelSize;
sFEM.polarTheta = (0:dt:(360-dt)) * (pi/180);


% Elliptic coordinate system
coefs = sFEM.coefsInit;
xa = sFEM.xa - coefs(1);
ya = sFEM.ya - coefs(2);
% Correction factors
if abs(coefs(4)) > 0
    p0 = -atan((1-coefs(3)+sqrt((coefs(3)-1)^2+coefs(4)^2))/coefs(4));
else
    p0 = 0;
end
a0 = sqrt(2*(1+coefs(3)+sqrt((coefs(3)-1)^2+coefs(4)^2)) ...
    / (4*coefs(3)-coefs(4)^2));
b0 = sqrt(2*(1+coefs(3)-sqrt((coefs(3)-1)^2+coefs(4)^2)) ...
    /(4*coefs(3)-coefs(4)^2));
ratio = b0 / a0;
m = [ratio*cos(p0)^2 + sin(p0)^2 ...
    -cos(p0)*sin(p0) + ratio*cos(p0)*sin(p0);
    -cos(p0)*sin(p0) + ratio*cos(p0)*sin(p0) ...
    cos(p0)^2 + ratio*sin(p0)^2];
ta = atan2(m(2,1)*xa + m(2,2)*ya, m(1,1)*xa + m(1,2)*ya);
ra = sqrt(xa.^2 + coefs(3)*ya.^2 + coefs(4)*xa.*ya) * b0;


% Resamping coordinates
Nout = [length(sFEM.polarRadius) length(sFEM.polarTheta)];
rInd = round((ra - sFEM.polarRadius(1)) / dr) + 1;
tInd = mod(round((ta - sFEM.polarTheta(1)) / (dt*pi/180)), Nout(2)) + 1;
sub = rInd <= Nout(1) & rInd >= 1;
rtIndsSub = [rInd(sub) tInd(sub)];

% KDE kernels
s = rSigma / dr;
v = ((-ceil(4*s)):ceil(4*s))';
kr = exp(-v.^2/(2*s^2));
s = tSigma / dt;
v = ((-ceil(4*s)):ceil(4*s));
kt = exp(-v.^2/(2*s^2));

% Normalization, mask and mean CBED
sFEM.polarNorm = accumarray(rtIndsSub,ones(sum(sub(:)),1),Nout);
sFEM.polarMask = accumarray(rtIndsSub,sFEM.mask_full(sub),Nout);
sFEM.polarCBEDmean = accumarray(rtIndsSub,sFEM.CBEDmean(sub),Nout);
% % KDE
kNorm = 1./ convolve2(conv2(ones(Nout),kr,'same'),kt,'wrap');
sFEM.polarNorm = kNorm .* convolve2(conv2( ...
    sFEM.polarNorm,kr,'same'),kt,'wrap');
sFEM.polarMask = kNorm .* convolve2(conv2( ...
    sFEM.polarMask,kr,'same'),kt,'wrap');
sFEM.polarCBEDmean = kNorm .* convolve2(conv2( ...
    sFEM.polarCBEDmean,kr,'same'),kt,'wrap');
% Normalize
sFEM.polarNormMask = sFEM.polarNorm > 0.5; % min intensity
sFEM.polarMask(sFEM.polarNormMask) = ...
    sFEM.polarMask(sFEM.polarNormMask) ...
    ./ sFEM.polarNorm(sFEM.polarNormMask);
sFEM.polarCBEDmean(sFEM.polarNormMask) = ...
    sFEM.polarCBEDmean(sFEM.polarNormMask) ...
    ./ sFEM.polarNorm(sFEM.polarNormMask);
sFEM.polarMask(~sFEM.polarNormMask) = 0;
sFEM.polarCBEDmean(~sFEM.polarCBEDmean) = 0;

% Loop through all images and compute polar transform
sFEM.polarAll = zeros(Nout(1),Nout(2),sFEM.stackSize(3));
for a0 = 1:sFEM.stackSize(3)
    CBED = sFEM.stack3D(:,:,a0);
    polarCBED = accumarray(rtIndsSub,CBED(sub),Nout);
    polarCBED = kNorm .* convolve2(conv2( ...
        polarCBED,kr,'same'),kt,'wrap');
    polarCBED(sFEM.polarNormMask) = ...
        polarCBED(sFEM.polarNormMask) ...
        ./ sFEM.polarNorm(sFEM.polarNormMask);
    polarCBED(~sFEM.polarNormMask) = 0;
    sFEM.polarAll(:,:,a0) = polarCBED;
    
    comp = a0 / sFEM.stackSize(3);
    progressbar(comp,2);
end

sFEM.p0 = p0;

% figure(11)
% clf
% % imagesc(sFEM.polarNorm)
% % imagesc(sFEM.polarMask)
% imagesc(sFEM.polarCBEDmean)
% colorbar



% % Full stack
% % if nargin > 1
% sFEM.polarCBEDvar = zeros(Nout);
% polarCBED = zeros(Nout);
% CBED = zeros(sFEM.stackSize(1:2));
% 
% for a0 = 1:size(stack3D,3)
%     CBED(:) = stack3D(:,:,a0);
%     polarCBED(:) = accumarray(rtIndsSub,CBED(sub),Nout);
%     
%     polarCBED = kNorm .* convolve2(conv2( ...
%         polarCBED,kr,'same'),kt,'wrap') ...
%         ./ sFEM.polarNorm;
%     
%     sFEM.polarCBEDvar(:) = sFEM.polarCBEDvar ...
%         + (polarCBED - sFEM.polarCBEDmean).^2;
% end
% sFEM.polarCBEDvar(:) = sFEM.polarCBEDvar / size(stack3D,3);
% 
% % Plot the mean variance over theta
% yMean = mean(sFEM.polarCBEDmean,2);
% yVar = mean(sFEM.polarCBEDvar,2);
% yNorm = mean(sFEM.polarMask,2);
% sub = yNorm > 1e-1;
% yMean(sub) = yMean(sub) ./ yNorm(sub);
% yVar(sub) = yVar(sub) ./ yNorm(sub);
% yVar(~sub) = 0;
% 
% yy = yVar;
% if flagNormPlot == true
%     sub = yMean > 0;
%     yy(sub) = yVar(sub) ./ yMean(sub).^2;
% end
% if flagPlot == true
%     figure(12)
%     clf
%     plot(sFEM.polarRadius,yy,'linewidth',2,'color','r')
%     if flagNormPlot == true
%         ylim(yRangeNorm)
%     else
%         ylim(yRange)
%     end
%     xlim(sFEM.polarRadius([1 end]))
%     xlabel('Scattering Vector [1/nm]')
% end
% % end
% 
% % Output results
% sFEM.yVar = yVar;
% sFEM.yMean = yMean;
% sFEM.yy = yy;
% 
% 
% if flagPlot == true
%     % plotting
%     figure(11)
%     clf
%     % Ip = sFEM.polarCBEDmean;
%     % Ip = kNorm;
%     if size(stack3D,3) > 1
%         Ip = sFEM.polarCBEDvar;
%     else
%         Ip = max(sFEM.polarCBEDmean,0).^0.5;    
%     end
%     imagesc(Ip,...
%         'xdata',sFEM.polarTheta*180/pi,...
%         'ydata',sFEM.polarRadius)
%     colormap(gray(256))
%     colorbar
% %     caxis(cRange)
% end
% % % Ip = (abs(ra - 244) < 2) .* ra;
% % % imagesc(Ip)
% % % axis equal off
% % % colormap(jet(256))
% % 
% % 
% % % % plotting
% % % figure(11)
% % % clf
% % % % Ip1 = (abs(ta - pi/4) < 0.04) ...
% % % %     +(abs(ta - pi*3/4) < 0.04) ...
% % % %     +(abs(ta + pi/4) < 0.04) ...
% % % %     +(abs(ta + pi*3/4) < 0.04);
% % % % Ip2 = (abs(ta - 0) < 0.04) ...
% % % %     +(abs(ta - pi/2) < 0.04) ...
% % % %     +(abs(ta + pi/2) < 0.04) ...
% % % %     +(abs(ta + pi) < 0.04) ...
% % % %     +(abs(ta - pi) < 0.04);
% % % % Ip = Ip2 - Ip1;
% % % % Ip = (abs(ra - 100) < 5) .* (mod(ta+0*pi/8,pi/4)-pi/8) * (180/pi);
% % % Ip = (abs(ra - 100) < 5) .* ra;
% % % imagesc(Ip)
% % % % imagesc(ra .* sin(ta))
% % % % imagesc(ra)
% % % axis equal off
% % % colorbar
% % % colormap(jet(256))
% % 
% % radiusReference = 130;
% % phi = 16 * pi/180;  % coordinate transform rotation value
% % yR = [-1 1]*0.4;
% % fsize = 18;
% % strainRef = [0 0 0] / 100;
% % 
% % % Rotate the coordinates, compute strains, display components
% % 
% % % Get fitted coefficients
% % %coefsAll = sStrain.coefsAll;
% % coefsAll = sFEM.fitCoefs;
% % r = coefsAll(:,10);
% % B = coefsAll(:,4);
% % C = coefsAll(:,3);
% % 
% % % Rotate coordinates
% % c2p = cos(2*phi);
% % s2p = sin(2*phi);
% % Arot = 0.5*(1 + C + (1-C)*c2p - B*s2p);
% % Brot = B*c2p + (1-C)*s2p;
% % Crot = 0.5*(1 + C + (C-1)*c2p + B*s2p);
% % 
% % 
% % % Scale to reference radius, i.e. normal ellipse form
% % a = Arot .* (radiusReference./r).^2;
% % b = 0.5 * Brot .* (radiusReference./r).^2;
% % c = Crot .* (radiusReference./r).^2;
% % 
% % 
% % % Convert to strains: minus sign out front = reciprocal space strain
% % exx = -0.5 * (a - 1);
% % eyy = -0.5 * (c - 1);
% % exy = -0.5 * b;
% % 
% % % subtract reference strains
% % exx = exx - strainRef(1);
% % eyy = eyy - strainRef(2);
% % exy = exy - strainRef(3);
% % 
% % 
% % % exxRMS = exx - mean(exx);
% % exxRMS = exx;
% % exxRMS = sqrt(mean(exxRMS.^2));
% % % eyyRMS = eyy - mean(eyy);
% % eyyRMS = eyy;
% % eyyRMS = sqrt(mean(eyyRMS.^2));
% % % exyRMS = exy - mean(exy);
% % exyRMS = exy;
% % exyRMS = sqrt(mean(exyRMS.^2));
% % 
% % sFEM.strains = [exxRMS eyyRMS exyRMS];


toc
end