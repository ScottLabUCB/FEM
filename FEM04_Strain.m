function [sFEM] = FEM04_Strain(sFEM)

% Colin Ophus and Ellis Kennedy - Sept 2019
% FEM analysis 04 - Strain fitting from elliptic coefficients; uses
% intrinsic ellipse equation to calculate dilation/expansion and deviation
% of major and minor axes from Cartesian axes; gives xx, yy, and xy strain
% percentages; plots data ellipse withrefernce fitting circle

%Inputs: 1 - sFEM from FEM03_Variance


% Reference ellipse
radiusReference = sFEM.fitCoefs(8);
phiReference = 0 * pi / 180;


% Get coefs
B = sFEM.fitCoefs(4);
C = sFEM.fitCoefs(3);
R = sFEM.fitCoefs(8);
R0 = R;

% Measure strain in distorted ellipse
if B == 0
    m11 = (R/radiusReference);
    m22 = (R/radiusReference) / sqrt(C);
    m12 = 0;
else
    tt = 4*C-B^2;
    den = sqrt((B^2+(C-1)^2)*tt/(1+C-sqrt(tt)));
    m11 = (R/radiusReference) * (2*C + sqrt(tt)) / den;
    m22 = (R/radiusReference) * (2 + sqrt(tt)) / den;
    m12 = -(R/radiusReference) * B / den;
end
mMeas = [m11 m12; m12 m22];
% Invert
mMeas = inv(mMeas);
% rotate into reference angle direction - Tensor rotation!
m = [cos(phiReference) sin(phiReference);
    -sin(phiReference) cos(phiReference)];
mMeas = m * mMeas * m';
% Extract strains
exxMeas = mMeas(1,1) - 1;
eyyMeas = mMeas(2,2) - 1;
exyMeas = (mMeas(1,2) + mMeas(2,1));
if nargout == 0
    disp(['exx = ' sprintf('%.03f',exxMeas*100) '%, ' ...
        'exy = ' sprintf('%.03f',exyMeas*100) '%, ' ...
        'eyy = ' sprintf('%.03f',eyyMeas*100) '%'])
else
    sFEM.strainMeas = [exxMeas*100 eyyMeas*100 exyMeas*100];
end


dt = 0.1;
t = (0:dt:(360-dt))' * pi/180;
xy0 = [cos(t) sin(t)] * R0;
m = [1+exxMeas exyMeas/2;
    exyMeas/2 1+eyyMeas];
m = inv(m);  % reciprocal space!
xy = xy0 * m;

figure(11)
clf
hold on
scatter(xy0(:,2),xy0(:,1),'k.')
scatter(xy(:,2),xy(:,1),'r.')
hold off
axis equal
legend('reference circle', 'fitted signal')
legend('boxoff')
set(gca,'ydir','reverse')


end