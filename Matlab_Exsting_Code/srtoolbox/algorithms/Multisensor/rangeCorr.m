function [cMult, cAdd, inlier] = rangeCorr(Iref, Itemp, method, varargin)
% Range data correction for range super-resolution.

    if nargin < 3
        method = 'linear';
    end
    
    % Get pixel values from reference image and template image.
    % We consider only non-zero pixel values within the images to suppress
    % "black borders" after geometric alignment of both images.
    xref = Iref( (Iref ~= 0) & (Itemp ~= 0) );
    xtemp = Itemp( (Iref ~= 0) & (Itemp ~= 0) );
    x(1,:) = xref;
    x(2,:) = xtemp;
    
    if strcmp(method, 'linear')    
        % Fit affine photometric model to the given pixel values by linear
        % approach.
        p = fitAffineModel(x); 
    elseif strcmp(method, 'msac')
        % Use MSAC to fit affine photometric model.
        if nargin < 4
            sigma = 0.025;
        else
            sigma = varargin{1}; 
        end
        if nargin < 5
            maxTrails = 500;
        else
            maxTrails = varargin{2};
        end
        [p, inlier] = msac(x, @fitAffineModel, @evalAffineModel, @checkDegenerationAffineModel, 2, 1.96*sigma, maxTrails, [1 0]);
        p = fitAffineModel(x(:, inlier));
    else      
        error('Unknown method %s', method);       
    end
    
    % Get multiplicative and additive part of the affine model.
    cMult = p(1);
    cAdd = p(2);
        
function p = fitAffineModel(x)

    p = polyfit(x(1,:), x(2,:), 1);

function [dist, inliers] = evalAffineModel(p, x, thresh)

    dist = (p(1) * x(1,:) + p(2)) - x(2,:);
    inliers = find(dist.^2 < thresh.^2);
    
function deg = checkDegenerationAffineModel(x)

    if abs(x(1,1) - x(1,2)) < eps
        deg = true;
        return;
    end
    
    p = fitAffineModel(x);
    if p(1) < 0
        % Degenerate solution: found negative slope
        deg = true;
        return;
    end
        
    deg = false;