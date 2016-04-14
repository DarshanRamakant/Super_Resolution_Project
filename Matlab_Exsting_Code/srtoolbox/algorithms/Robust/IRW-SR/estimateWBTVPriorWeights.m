function [weights, scaleParameter] = estimateWBTVPriorWeights(z, p, weights, scaleParameter)

    if nargin < 2
        % Use default value p = 0.5 for the sparsity parameter.
        p = 0.5;
    end
    
    if nargin < 4
        % Adaptive estimation of the scale parameter.
        scaleParameter = getAdaptiveScaleParameter(z, weights);
    end
    
    % Estimation of the weights based on pre-selected scale parameter.
    c = 2;
    scaleParameter = c * scaleParameter;
    weights = (p * scaleParameter^(1-p)) ./ (abs(z).^(1-p));
    weights(abs(z) <= scaleParameter) = 1;
    
function scaleParameter = getAdaptiveScaleParameter(z, weights)

    scaleParameter = weightedMedian( abs(z - weightedMedian(z, weights)), weights );



