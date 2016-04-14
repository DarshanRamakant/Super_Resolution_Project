function [SR_input, SR_guide] = guidedsr(I_input, model_input, I_guide, model_guide, weightFcn, gf_radius, gf_eps, gf_lambda, varargin)
    
    % Set default parameters for guided filtering if no user-defined
    % parameters are available.
    if nargin < 6
        gf_radius = 1;
    end
    if nargin < 7
        gf_eps = 1e-4;
    end
    if nargin < 8
        gf_lambda = 0.5;
    end
    
    if nargin > 8
        solverParams = varargin{end};
    else
        % Use default solver parameters
        solverParams = SRSolverParams;
    end
     
    % Assemble equation system to be solved consisting of low-resolution
    % observations and the system matrix.
    lrDim = size(I_input(:,:,1));
    gdDim = size(I_guide(:,:,1));
    if solverParams.verbose
        disp('Assemble equation system (system matrix)...');
    end
    [W_input, I_input] = composeSREquationSystem(I_input, model_input);
    [W_guide, I_guide] = composeSREquationSystem(I_guide, model_guide);
    
    if isempty(model_input.SR)
        % Compute initial guess for the super-resolved image. We use the
        % "average image" estimated from the system matrix and photometric
        % parameters (if available).
        if solverParams.verbose
            disp('Compute average image used as initial guess...');
        end
        SR_input = getInitialSRImage(W_input, I_input, model_input.photometricParams);
    else
        % Use initial guess provided by the user and reshape 2-D image into
        % parameter vector.
        SR_input = imageToVector(model_input.SR);
    end
        
    if isempty(model_guide.SR)
        % Compute initial guess for the super-resolved image. We use the
        % "average image" estimated from the system matrix and photometric
        % parameters (if available).
        if solverParams.verbose
            disp('Compute average image used as initial guess...');
        end
        SR_guide = getInitialSRImage(W_guide, I_guide, model_guide.photometricParams);
    else
        % Use initial guess provided by the user and reshape 2-D image into
        % parameter vector.
        SR_guide = imageToVector(model_guide.SR);
    end
     
    scgOptions = setupSCGOptions(solverParams);
    if solverParams.verbose
        disp('Minimize objective function...');
    end
    
    conf_guide = imageToVector(model_guide.confidence);
    if isempty(conf_guide)
        conf_guide = ones(size(I_guide));
    end
    conf_input = imageToVector(model_input.confidence);
    if isempty(conf_input)
        conf_input = ones(size(I_input));
    end
    
    % Perform IRLS iterations.
    for iter = 1 : solverParams.maxIrlsIter
        
        % Step 1: Determine confidence weights
        % Get residual error for the guidance images.
        r_guide = getResidual(SR_guide, I_guide, W_guide, model_guide.photometricParams);
        weights_guide = weightFcn(r_guide);
        % Combine IRLS weights with user-defined confidence map.
        model_guide.confidence = conf_guide .* weights_guide;
        % Get residual error for the input images.
        r_input = getResidual(SR_input, I_input, W_input, model_input.photometricParams);
        weights_input = weightFcn(r_input);
        % Combine IRLS weights with user-defined confidence map.
        model_input.confidence = conf_input .* weights_input;
                 
        % Step 2: Super-resolution for guidance images
        SR_guide_prev = SR_guide;
        SR_guide = scg(@mapfunc, SR_guide', scgOptions, @mapfunc_grad, model_guide, I_guide, W_guide, W_guide');      
        SR_guide = SR_guide';
        
        % Step 3: Updated guided filter coefficients
        I = vectorToImage(SR_guide, model_guide.magFactor * gdDim);
        p = vectorToImage(SR_input, model_input.magFactor * lrDim); 
        [A, b] = guidedfilter(I, p, gf_radius, gf_eps); 
        % Reorganize filter coefficents for fast matrix vector computation.
        A = spdiags(imageToVector(A), 0, numel(A), numel(A));
        b = imageToVector(b);
    
        % Step 4: Super-resolution for input imates
        SR_input_prev = SR_input;
        SR_input = scg(@mapfuncGuided, SR_input', scgOptions, @mapfuncGuided_grad, model_input, I_input, W_input, W_input', A, b, SR_guide, gf_lambda); 
        SR_input = SR_input';
        
        % Test for convergence.
        if max(abs(SR_guide_prev - SR_guide)) < solverParams.tolX && max(abs(SR_input_prev - SR_input)) < solverParams.tolX
            % Image reconstruction converged.
            break;
        end
            
    end
    
    % Reshape parameter vector to a 2D image.
    SR_input = vectorToImage(SR_input, model_input.magFactor * lrDim);
    SR_guide = vectorToImage(SR_guide, model_guide.magFactor * gdDim);
    

function scgOptions = setupSCGOptions(solverParams)

    scgOptions = zeros(1,18); 
    scgOptions(1) = solverParams.verbose;        
    scgOptions(2) = solverParams.tolX;   
    scgOptions(3) = solverParams.tolF;
    scgOptions(9) = solverParams.gradCheck;
    scgOptions(10) = solverParams.maxFunEvals;     
    scgOptions(14) = solverParams.maxIter;

function f = mapfuncGuided(SR, model, LR, W, Wt, A, b, SR_guide, gf_lambda)
    
    if ~iscolumn(SR)
        % Reshape to column vector. 
        SR = SR';
    end
    
    % Evaluate the data fidelity term.
    dataTerm = mapDataTerm(SR, model, LR, W);
    
    % Evaluate image prior for regularization the super-resolved estimate.
    prior = model.imagePrior.function(SR, model.imagePrior.parameters{1:end});
    
    % Evaluate data fidelity for guided filter
    dataTerm_gf = sum((SR - A*SR_guide - b).^2);
       
    f = dataTerm + gf_lambda * dataTerm_gf + model.imagePrior.weight * prior;
   
function grad = mapfuncGuided_grad(SR, model, LR, W, Wt, A, b, SR_guide, gf_lambda)
    
    if ~iscolumn(SR)
        % Reshape to column vector. 
        SR = SR';
    end
    
    % Calculate gradient of the data fidelity term w.r.t. the
    % super-resolved image.
    dataTerm_grad = mapDataTerm_gradImage(SR, model, LR, W, Wt);
    
    dataTerm_gf_grad = 2*(SR - A*SR_guide - b);
    
    % Calculate gradient of the image prior w.r.t. the super-resolved image.
    prior_grad = model.imagePrior.gradient(SR, model.imagePrior.parameters{1:end});
    
    % Sum up to total gradient
    grad = dataTerm_grad + gf_lambda * dataTerm_gf_grad + model.imagePrior.weight * prior_grad;
    grad = grad';
    
    
function f = mapfunc(SR, model, LR, W, ~)
    
    if ~iscolumn(SR)
        % Reshape to column vector. 
        SR = SR';
    end
    
    % Evaluate the data fidelity term.
    dataTerm = mapDataTerm(SR, model, LR, W);
    
    % Evaluate image prior for regularization the super-resolved estimate.
    prior = model.imagePrior.function(SR, model.imagePrior.parameters{1:end});
       
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
