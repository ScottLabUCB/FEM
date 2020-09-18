function [sFEM] = FEM42(sFEM)

% Colin Ophus - 2019 Sept
% 42 - normalize all CBED images inside ring mask


sFEM.ringInt = zeros(sFEM.stackSize(3),1);
for a0 = 1:sFEM.stackSize(3)
    I = sFEM.stack3D(:,:,a0);
    sFEM.ringInt(a0) = mean(I(sFEM.mask_ring)); %mean
end

sFEM.scaleInt = mean(sFEM.ringInt) ./ sFEM.ringInt;  %mean
for a0 = 1:sFEM.stackSize(3)
    sFEM.stack3D(:,:,a0) = ...
        sFEM.stack3D(:,:,a0) * sFEM.scaleInt(a0);
end
sFEM.CBEDmean = median(sFEM.stack3D,3);  %mean or median


% % plotting for testing
% Ip = reshape(sFEM.ringInt,[12 12*5]);
% figure(11)
% clf
% imagesc(Ip)
% axis equal off
% colorbar


end

