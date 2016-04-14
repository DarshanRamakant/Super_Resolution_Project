function opts = getReweightedOptimizationParams

    % Maximum number of majorization-minimization (MM) iterations used for
    % iteratively re-weighted minimization (should be between 5 and 15).
    opts.maxMMIter = 10;
    
    % Maximum number of scaled conjugate gradients (SCG) iterations in the
    % inner optimization loop used for image reconstruction (should be
    % between 3 and 10).
    opts.maxSCGIter = 5;
    
    % Maximum number of cross validation (CV) iterations for hyperparameter
    % selection (should be between 10 and 30).
    opts.maxCVIter = 20;
        
    % Fraction of the observations used for the training stage in
    % hyperparameter selection (should be between 0.90 and 0.95).
    opts.fractionCVTrainingObservations = 0.95;
    
    % Initial search range (lower and upper bound) for the regularization 
    % hyperparameter on a logarithmic scale (should be between 10-12 and
    % 10^0).
    opts.hyperparameterCVSearchRange = [-12 0];
    
    % Termination tolerance for max(x(t) - x(t-1)) for the estimate of the
    % high-resolution image over the MM and SCG iterations (should be 
    % approx. 1e-3).
    opts.terminationTol = 1e-3;
    
    % Sparsity parameter used for weighted bilateral total variation prior
    % (should be approx. 0.5 and can be selected in [0; 1]).
    opts.sparsityParameter = 0.5;
    
    % Number of levels used for coarse-to-fine optimization. If empty, the
    % number of levels is set automatically (should be between 2 and 5).
    opts.numCoarseToFineLevels = [];
    
    