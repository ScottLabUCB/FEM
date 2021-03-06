function [sFEM] = FEM43(sFEM,coefs)
tic

% Colin Ophus - May 2019
% 43 fit elliptic coordinate system


flagReportError = 0;
flagPlot = 1;
powerPlot = 0.5;
skipFit = [11 1];

% powerLawRing = 3.0;

if ~isfield(sFEM,'coefsInit')
    % Ellis experiments Ta / Si3N4
    
%     % Ta
%     sFEM.coefsInit = [ ...
%         255 254 0.96 0.07 ...
%         0 ...
%         930 94 ...
%         121 370 9 9.5];
%     
%      % SiN
%     sFEM.coefsInit = [ ...
%         236 250 0.96 0.07 ...
%         50 ...
%         1000 70 ...
%         121 100 30 30];
    
       %NbO
%     sFEM.coefsInit = [ ...
%         285 250 0.99 0.01 ...
%         30 ...
%         900 80 ...
%         130 100 20 20];
    
     %RT
%     sFEM.coefsInit = [ ...
%         270 260 1.00 0.00 ...
%         200 ...
%         10000 20 ...
%         150 200 15 15];
    
    % 300 SiN
    sFEM.coefsInit = [ ...
        260 270 1.00 0.00 ...
        200 ...
        10000 20 ...
        150 200 15 15 ...
        ];

end
if ~isfield(sFEM,'basis')
    sFEM.basis = [sFEM.xa(:) sFEM.ya(:)];
end
lb = [0 0 0.5 -0.5 ...
    0 0 0 0 0 0 0];
ub = [sFEM.stackSize(1:2) 2 0.5 ...
    Inf Inf Inf Inf Inf Inf Inf];

sFEM.fitOptions = optimset( ...
    'TolX',1e-12,...
    'TolFun',1e-12,...
    'MaxFunEvals',1e4,...
    'maxiter',1e4,...
    'display','off');
sFEM.fitFun = @(c,x) c(5) ...
    ...
    + c(6)*exp((-1/2/c(7)^2) ...
    *((x(:,1) - c(1)).^2 ...
    + c(3)*(x(:,2) - c(2)).^2 ...
    + c(4)*(x(:,1) - c(1)).*(x(:,2) - c(2)))) ...
    ...
    + c(9)*exp((-1/2/c(10)^2)* ...
      abs(c(8) - sqrt( ...
      (x(:,1) - c(1)).^2 ...
    + c(3)*(x(:,2) - c(2)).^2 ...
    + c(4)*(x(:,1) - c(1)).*(x(:,2) - c(2)))).^2) ...
    .* (c(8)^2 ...
    > ((x(:,1) - c(1)).^2 ...
    + c(3)*(x(:,2) - c(2)).^2 ...
    + c(4)*(x(:,1) - c(1)).*(x(:,2) - c(2)))) ...
    ...
    + c(9)*exp((-1/2/c(11)^2)* ...
      abs(c(8) - sqrt( ...
      (x(:,1) - c(1)).^2 ...
    + c(3)*(x(:,2) - c(2)).^2 ...
    + c(4)*(x(:,1) - c(1)).*(x(:,2) - c(2)))).^2) ...
    .* (c(8)^2 ...
    < ((x(:,1) - c(1)).^2 ...
    + c(3)*(x(:,2) - c(2)).^2 ...
    + c(4)*(x(:,1) - c(1)).*(x(:,2) - c(2))));


% If output is specified, perform fitting
if nargout > 0 && nargin < 2
    inds2D = sub2ind(sFEM.stackSize(1:2),sFEM.xa,sFEM.ya);
    
    % First fit
    for a0 = 1:length(skipFit)
        maskFit = mod(inds2D,skipFit(a0)) == 0;
        
        sFEM.coefsInit = lsqcurvefit(sFEM.fitFun, ...
            sFEM.coefsInit,...
            sFEM.basis(sFEM.mask_ring(:) & maskFit(:),:),...
            double(sFEM.CBEDmean(sFEM.mask_ring(:) & maskFit(:))),...
            lb,ub,sFEM.fitOptions);
    end
end
if nargin == 2
    sFEM.coefsInit = coefs;
end




% Plotting
if flagPlot == true
    Imeas = max(sFEM.CBEDmean,0);
    Ifit = reshape(sFEM.fitFun(sFEM.coefsInit,sFEM.basis), ...
        sFEM.stackSize(1:2)); 
    
    Ip1 = Imeas.^powerPlot;
    Ip2 = Ifit .^powerPlot;
    t = atan2(sFEM.ya - sFEM.coefsInit(2),...
        sFEM.xa - sFEM.coefsInit(1));
    w = double(mod(t,pi/8) < pi/16);
    w(~sFEM.mask_ring & (w == 1)) = 0;
    Ip = Ip1 .* w ...
        + Ip2 .* (1-w);
    intRange = [min(Ip1(sFEM.mask_ring)) ...
        max(Ip1(sFEM.mask_ring))];
    
    figure(11)
    clf
    imagesc(Ip)
    axis equal off
    colormap(jet(256))
    set(gca,'position',[0 0 1 1])
    caxis(intRange)
    drawnow;
    % %     cR = [min(Ip1(sFEM.mask)) max(Ip1(sFEM.mask))];
    % %     caxis([cR(1) 0.6*cR(1) + 0.4*cR(2)])
    %     caxis([10 40])
    
    
    if flagReportError == true
        fitError = mean(abs( ...
            Imeas(sFEM.mask_ring) - Ifit(sFEM.mask_ring))) ...
            / mean(abs(Imeas(sFEM.mask_ring)));
        disp(['Mean abs. diff. = ' ...
            sprintf('%.04f',fitError*100) '%'])
        
    end
end



toc
end