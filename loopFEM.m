function [sFEM] = loopFEM(stack3D,sROIfull,coefs)

sFEM = FEM01(stack3D,sROIfull);

% If user provides coefs, use as initial guess
if nargin > 3
    sFEM.fitCoefs = coefs;
end

%FEM02(sFEM);
sFEM = FEM02(sFEM);
sFEM = FEM03(sFEM,stack3D);

sFEM = FEM04(sFEM);

end