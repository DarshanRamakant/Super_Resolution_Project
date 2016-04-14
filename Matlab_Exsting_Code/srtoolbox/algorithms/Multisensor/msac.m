function [bestModel, bestInliers] = msac(x, fittingFun, distFun, degFun, minSamples, thresh, maxTrials, baselineModel, varargin)
% MSAC - M-estimator sample consensus estimation
%
% Usage:
%   [bestModel, bestInliers] = msac(x, fittingFun, distFun, minSamples, ...
%                                   thresh, maxTrials) 
%
% MSAC estimates a model free of outliers from given input data X according
% to a user-defined fitting and distance function. In contrast to RANSAC
% the cost function optimized in MSAC takes the error for inliers into 
% accoint (instead of assign constant cost to each inlier).
% See P. Torr and A. Zisserman, MLESAC: A New Robust Estimator with 
% Application to Estimating Image Geometry, 2000 for details.
%
% This implementation is in the style of the RANSAC code of Peter Kovesis 
% Matlab toolbox see http://www.csse.uwa.edu.au/~pk/research/MatlabFns/.
%
% Parameters:
%   x:          A n x m matrix with m samples of dimension n
%   fittingFun: Function handle to determine a model M = fittingFun(x)
%   distFun:    Distance function and inlier classification according to
%               [dist, inliers] = distFun(model, x, thresh)
%   minSamples: The minimum number of samples required to fit a model
%   thresh:     Threshold to classify for inliers/outliers
%   maxTrails:  The maximum number of trails used for estimation (default:
%               1000)
%
% Return:
%   bestModel:      The estimated model determined having minimum costs in
%                   MSAC algorithm
%   bestInliers:    The array of indices that are the inliers for estimated 
%                   model

    if nargin < 7
        maxTrials = 1000;
    end
    numTrials = 0;  
    N = maxTrials;
    numSamples = size(x, 2);
    
    if nargin < 8
        bestModel = [];
        bestCost = Inf;
    else
        % Setup baseline model and initialize optimal values.
        bestModel = baselineModel;
        [dist, inliers] = distFun(baselineModel, x, thresh, varargin{1:end});
        bestCost = evalCosts(dist, inliers, thresh);
    end
    
    bestInliers = 1:numSamples;
    
    while numTrials < N
        
        % Select 'minSamples' points from the given data 'x'.
        sampleIdx = randsample(numSamples, minSamples);
        
        if degFun(x(:,sampleIdx), varargin{1:end})
            % Found degenerate set. Proceed with next iterations.
            numTrials = numTrials + 1;
            continue;
        end
        
        % Fit a model to the selected subset of samples
        model = fittingFun( x(:,sampleIdx), varargin{1:end} );
          
        % Calculate distance for each point to the model and classify for
        % inliers according to the given threshold.
        [dist, inliers] = distFun(model, x, thresh, varargin{1:end});
        
        % Evalute MSAC cost function
        totalCost = evalCosts(dist, inliers, thresh);
        
        if totalCost < bestCost    
            % We found a model with minimum cost. Save inliers and the
            % estimated model.
            bestCost = totalCost;
            bestModel = model;
            bestInliers = inliers;
        end
        
        numTrials = numTrials + 1;
        if numTrials > maxTrials
            % User-defined maximum number of trails reached.
            break;
        end
        
    end
    
 function costs = evalCosts(dist, inliers, thresh)    
 
     costs = zeros(size(dist));
     costs(1:end) = thresh;           % Costs for outliers
     costs(inliers) = dist(inliers);  % Costs for inliers
     costs = sum(costs.^2);