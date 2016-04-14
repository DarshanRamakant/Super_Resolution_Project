function SR = irlssr(LRImages, model, weightFcn, varargin)
    
    if nargin > 3
        solverParams = varargin{end};
    else
        % Use default solver parameters
        solverParams = SRSolverParams;
    end
     
    % Assemble equation system to be solved consisting of low-resolution
    % observations and the system matrix.
    lrDim = size(LRImages(:,:,1));
    if solverParams.verbose
        disp('Assemble equation system (system matrix)...');
    end
    [W, LRImages] = composeSREquationSystem(LRImages, model);
        
    if isempty(model.SR)
        % Compute initial guess for the super-resolved image. We use the
        % "average image" estimated from the system matrix and photometric
        % parameters (if available).
        if solverParams.verbose
            disp('Compute average image used as initial guess...');
        end
        SR = getInitialSRImage(W, LRImages, model.photometricParams);
    else
        % Use initial guess provided by the user and reshape 2-D image into
        % parameter vector.
        SR = imageToVector(model.SR);
    end
     
    scgOptions = setupSCGOptions(solverParams);
    if solverParams.verbose
        disp('Minimize objective function...');
    end
    
    iter = 0;
    SR_prev = Inf(size(SR));
    confWeights = imageToVector(model.confidence);
    if isempty(confWeights)
        confWeights = ones(size(LRImages));
    end
    while (iter < solverParams.maxIrlsIter) && (max(abs(SR_prev - SR)) > solverParams.tolX)
        
        % Get residual error
        r = getResidual(SR, LRImages, W, model.photometricParams);
        % Evaluate weight function
        weights = weightFcn(r);
        % Compose total weights
        model.confidence = weights .* confWeights;
        
        % Obtain refined solution for the HR image
        SR_prev = SR;
        SR = scg(@mapfunc, SR', scgOptions, @mapfunc_grad, model, LRImages, W, W');
        SR = SR';
                        
        iter = iter + 1;
    
    end
        
    % Reshape parameter vector to a 2D image.
    SR = vectorToImage(SR, model.magFactor * lrDim);
    

function scgOptions = setupSCGOptions(solverParams)

    scgOptions = zeros(1,18); 
    scgOptions(1) = solverParams.verbose;        
    scgOptions(2) = solverParams.tolX;   
    scgOptions(3) = solverParams.tolF;
    scgOptions(9) = solverParams.gradCheck;
    scgOptions(10) = solverParams.maxFunEvals;     
    scgOptions(14) = solverParams.maxIter;

function f = mapfunc(SR, model, LR, W, ~)
    
    if ~iscolumn(SR)
        % Reshape to column vector. 
        SR = SR';
    end
    
    % Evaluate the data fidelity term.
    dataTerm = mapDataTerm(SR, model, LR, W);
    
    % Evaluate image prior for regularization the super-resolved estimate.
    prior = model.imagePrior.function(SR, model.imagePrior.parameters{1:end});
    
    % Calculate objective function.
    f = dataTerm + model.imagePrior.weight * prior;
                
function grad = mapfunc_grad(SR, model, LR, W, Wt)
    
    if ~iscolumn(SR)
        % Reshape to column vector. 
        SR = SR';
    end
    
    % Calculate gradient of the data fidelity term w.r.t. the
    % super-resolved image.
    dataTerm_grad = mapDataTerm_gradImage(SR, model, LR, W, Wt);
    
    % Calculate gradient of the regularization term w.r.t. the 
    % super-resolved image.
    prior_grad = model.imagePrior.gradient(SR, model.imagePrior.parameters{1:end});
    
    % Sum up to total gradient
    grad = dataTerm_grad + model.imagePrior.weight * prior_grad;
    grad = grad';