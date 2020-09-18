function [sFEM] = FEM46(sFEM, phiReference)

% Colin Ophus and Ellis Kennedy - 2019 Sept
% 46 - strain

%%RT 0 tilt mean radius: 139.41; 300 0 tilt mean radius: 140.17
radiusReference = 139.41;  %was 122.725
if nargin == 1
    phiReference = 0* pi / 180; %%set to 0 if unknown
else
    phiReference = phiReference;
end

% Get coefs
B = sFEM.coefsInit(4);
C = sFEM.coefsInit(3);
R = sFEM.coefsInit(8);
% B = sFEM.fitCoefs(4);
% C = sFEM.fitCoefs(3);
% R = sFEM.fitCoefs(8);



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
exyMeas = (mMeas(1,2) + mMeas(2,1))/2;
if nargout == 0
    disp(['exx = ' sprintf('%.03f',exxMeas*100) '%, ' ...
        'exy = ' sprintf('%.03f',exyMeas*100) '%, ' ...
        'eyy = ' sprintf('%.03f',eyyMeas*100) '%'])
else
    sFEM.strainMeas = [exxMeas eyyMeas exyMeas];
    sFEM.mMeas = mMeas;
end


end

