function [sFEM] = FEM41(stack3D,sROI_full,sROI_ring)

% Colin Ophus - 2019 Sept
% FEM analysis 41 - initialize coordinte system, mask for beamstop, etc.

%stack3D = stack3D.cube;
sFEM.CBEDmean = mean(stack3D,3);
sFEM.stackSize = size(stack3D);
sFEM.stack3D = double(stack3D);
[sFEM.ya,sFEM.xa] = meshgrid(1:sFEM.stackSize(2),1:sFEM.stackSize(1));
%sFEM.scale = stack3D.scale(1);

% masking
sFEM.mask_ring = inpolygon(sFEM.xa,sFEM.ya,sROI_ring.p(:,1),sROI_ring.p(:,2));
sFEM.mask_full = inpolygon(sFEM.xa,sFEM.ya,sROI_full.p(:,1),sROI_full.p(:,2));


end