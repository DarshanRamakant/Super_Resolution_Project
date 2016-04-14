function [weights, scaleParameter] = estimateObservationConfidenceWeights(r, weights, scaleParameter)
          
    % Detect single frames without zero-mean residual as outliers.
    weightsVec = [];
    rVec = [];
    rMax = 0.02;
    for k = 1:length(r)   
        rVec = [rVec; r{k}];
        if abs( median(r{k}) ) < rMax
            weightsBias(k) = 1;
        else
            weightsBias(k) = 0;
        end   
        weightsVec = [weightsVec; weights{k} .* weightsBias(k)];
    end

    if nargin < 3
        % Use adaptive scale parameter estimation.
        scaleParameter = getAdaptiveScaleParameter(rVec(weightsVec > 0), weightsVec(weightsVec > 0));
    end

    for k = 1:length(r)

        % Estimate local confidence weights for current frame.
        c = 2;
        weightsLocal = 1 ./ abs(r{k});
        weightsLocal(abs(r{k}) < c*scaleParameter) = 1 / (c*scaleParameter);
        weightsLocal = c*scaleParameter * weightsLocal;

        % Assemble confidence weights from bias (frame-wise) weights and
        % local (pixel-wise) weights.
        weights{k} = weightsBias(k) .* weightsLocal;
        
    end
    
function scaleParameter = getAdaptiveScaleParameter (r, weights)

    scaleParameter = 1.4826 * weightedMedian( abs(r(weights > 0) - weightedMedian(r(weights > 0), weights(weights > 0))), weights(weights > 0) );

    