function [sFEM] = FEM01_MaskCoords(stack3D,sROIfull,sROIring)

% Colin Ophus - May 2019
% FEM analysis 01 - initialize coordinte system, applies full and ring...
% masks to remove beamstop.

% Inputs: 1 - image stack read in using dm4Reader; 2 - mask that excludes beamstop...
% and central bright spot from RealspaceLattice01; 3 - mask that excludes beamstop, central bright spot
% and everything past first ring from RealspaceLattice01

sFEM.CBEDmean = double(mean(stack3D,3));
sFEM.stackSize = size(stack3D);
[sFEM.ya,sFEM.xa] = meshgrid(1:sFEM.stackSize(2),1:sFEM.stackSize(1));

% masking
if nargin == 1
    sFEM.mask = true(sFEM.stackSize(1:2));
else
    sFEM.mask = inpolygon(sFEM.xa,sFEM.ya,sROIring.p(:,1),sROIring.p(:,2));
    sFEM.maskFull = inpolygon(sFEM.xa,sFEM.ya,sROIfull.p(:,1),sROIfull.p(:,2));
end



end