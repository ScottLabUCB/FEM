function [sFEM] = loopFEM04(stack3D,sROIfull,sROIring,coefs)

% Ellis Kennedy - April 2020
% FEM loop - Loops through FEM steps 41 - 46

sFEM = FEM41(stack3D,sROIfull,sROIring);

% If user provides coefs, use as initial guess
if nargin > 3
    sFEM.fitCoefs = coefs;
end

sFEM = FEM42(sFEM);
sFEM = FEM43(sFEM);
sFEM = FEM44(sFEM);
sFEM = FEM45(sFEM);
sFEM = FEM46(sFEM);

end