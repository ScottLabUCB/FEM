function [sFEM] = FEM45(sFEM)

% Colin Ophus - 2019 Sept
% variance

% Compute mean and variance
% polarCBEDmean = mean(sFEM.polarAll,3);
sFEM.polarCBEDvar = ...
    mean((sFEM.polarAll - sFEM.polarCBEDmean).^2,3);


% Compute radial stats
sFEM.radialNorm = sum(sFEM.polarMask,2);
sFEM.radialMean = sum(sFEM.polarCBEDmean,2);
sFEM.radialVar = sum(sFEM.polarCBEDvar,2);
% Normalzie
sub = sFEM.radialNorm > 1;
sFEM.radialMean(sub) = sFEM.radialMean(sub) ...
    ./ sFEM.radialNorm(sub);
sFEM.radialVar(sub) = sFEM.radialVar(sub) ...
    ./ sFEM.radialNorm(sub);
sFEM.radialMean(~sub) = 0;
sFEM.radialVar(~sub) = 0;
% Output
sFEM.radialVarNorm = sFEM.radialVar;
sFEM.radialVarNorm(sub) = ...
    sFEM.radialVarNorm(sub) ./ sFEM.radialMean(sub).^2;

% Ip = sFEM.polarCBEDvar ./ sFEM.polarCBEDmean.^2;

% plotting
figure(14)
clf
hold on
plot(sFEM.polarRadius,sFEM.radialVarNorm,...
    'linewidth',2,'color','r')
% imagesc(Ip .* sFEM.polarMask)
% imagesc([polarCBEDmean sFEM.polarCBEDmean])

xlim([0 sFEM.polarRadius(end)])
ylim([0 0.08])

end