function [sFEM] = FEM03_Variance(sFEM,stack3D)

% Colin Ophus - May 2019
% FEM analysis 03 - Polar coordinate transformation from elliptic
% coordinate system. Calculates variance and stores in sFEM struct.

%Inputs: 1 - SFEM from FEM02_FitFunction; 2 - dm4 stack after opening with
%dm4Reader

pixelSize = stack3D.scale(1); 
rSigma = 0.1;
tSigma = 2;

rMax = 220;
dr = 1; %1  
dt = 1; %5
sFEM.polarRadius = (0:dr:rMax)*pixelSize;
sFEM.polarTheta = (0:dt:(360-dt)) * (pi/180);

scale = 15 / 3 * 3;
yRangeNorm = [0 0.005] * scale;
yRange = [0 5e3] * scale;
cRange = [0 5e3] * scale;
flagPlot = true;
flagNormPlot = true;

% Elliptic coordinate system
coefs = sFEM.fitCoefs
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
% end
% [a0 b0]
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
sFEM.polarMask = accumarray(rtIndsSub,sFEM.maskFull(sub),Nout);
sFEM.polarCBEDmean = accumarray(rtIndsSub,sFEM.CBEDmean(sub),Nout);
% KDE
kNorm = 1./ convolve2(conv2(ones(Nout),kr,'same'),kt,'wrap');
sFEM.polarNorm = kNorm .* convolve2(conv2( ...
    sFEM.polarNorm,kr,'same'),kt,'wrap');
sFEM.polarMask = kNorm .* convolve2(conv2( ...
    sFEM.polarMask,kr,'same'),kt,'wrap') ...
    ./ sFEM.polarNorm;
sFEM.polarCBEDmean = kNorm .* convolve2(conv2( ...
    sFEM.polarCBEDmean,kr,'same'),kt,'wrap') ...
    ./ sFEM.polarNorm;

% Full stack
% if nargin > 1
sFEM.polarCBEDvar = zeros(Nout);
polarCBED = zeros(Nout);
CBED = zeros(sFEM.stackSize(1:2));

for a0 = 1:size(stack3D,3)
    CBED(:) = stack3D(:,:,a0);
    polarCBED(:) = accumarray(rtIndsSub,CBED(sub),Nout);
    
    polarCBED = kNorm .* convolve2(conv2( ...
        polarCBED,kr,'same'),kt,'wrap') ...
        ./ sFEM.polarNorm;
    
    sFEM.polarCBEDvar(:) = sFEM.polarCBEDvar ...
        + (polarCBED - sFEM.polarCBEDmean).^2;
end
sFEM.polarCBEDvar(:) = sFEM.polarCBEDvar / size(stack3D,3);

%%% for single CBEDmean
%      CBED = stack3D;
%      polarCBED = accumarray(rtIndsSub,CBED(sub),Nout);
%      
%      polarCBED = kNorm .* convolve2(conv2( ...
%          polarCBED,kr,'same'),kt,'wrap') ...
%          ./ sFEM.polarNorm;
%      
%      sFEM.polarCBEDvar = sFEM.polarCBEDvar ...
%          + (polarCBED - sFEM.polarCBEDmean).^2;



% Plot the mean variance over theta
yMean = mean(sFEM.polarCBEDmean,2);
yVar = mean(sFEM.polarCBEDvar,2);
yNorm = mean(sFEM.polarMask,2);
sub = yNorm > 1e-1;
yMean(sub) = yMean(sub) ./ yNorm(sub);
yVar(sub) = yVar(sub) ./ yNorm(sub);
yVar(~sub) = 0;

yy = yVar;
if flagNormPlot == true
    sub = yMean > 0;
    yy(sub) = yVar(sub) ./ yMean(sub).^2;
end
figure(12)
clf
plot(sFEM.polarRadius,yy,'linewidth',2,'color','r')
if flagNormPlot == true
    %ylim(yRangeNorm)
    ylim([0 0.5])
else
    %ylim(yRange)
    ylim([0 0.5])
end
xlim(sFEM.polarRadius([1 end]))
xlabel('Scattering Vector [1/nm]')
% end

% Output results
sFEM.yVar = yVar;
sFEM.yMean = yMean;
sFEM.yy = yy;

if flagPlot == true
    % plotting
    figure(11)
    clf
    % Ip = sFEM.polarCBEDmean;
    % Ip = kNorm;
    if size(stack3D,3) > 1
        Ip = sFEM.polarCBEDvar;
    else
        Ip = max(sFEM.polarCBEDmean,0).^0.5;    
    end
    imagesc(Ip,...
        'xdata',sFEM.polarTheta*180/pi,...
        'ydata',sFEM.polarRadius)
    colormap(gray(256))
    colorbar
%     caxis(cRange)
end

end