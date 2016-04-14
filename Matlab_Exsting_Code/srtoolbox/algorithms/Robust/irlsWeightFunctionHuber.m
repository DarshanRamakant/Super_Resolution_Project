function weights = irlsWeightFunctionHuber(r, thresh)
    
    if nargin < 2
        % Use adaptive threshold selection.
        thresh = 1.4826 * mad(r, 1);
    end
    
    % Weight is proportional to inverse of the residual error
    weights = 1 ./ abs(r);
    % Constant weight for small residual errors.
    weights(abs(r) < thresh) = 1 / thresh;
    weights(abs(r) > 5*thresh) = 0;
    
    % Normalization
    weights = thresh * weights;