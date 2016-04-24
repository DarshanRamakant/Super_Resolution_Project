function SR = mapsr(LRImages, model, varargin)
    % Use default solver parameters
    solverParams = SRSolverParams;
    
    % Assemble equation system to be solved consisting of low-resolution
    % observations and the system matrix.
    lrDim = size(LRImages(:,:,1));
    disp('Assemble equation system (system matrix)...');
   [W, LRImages] = composeSREquationSystem(LRImages, model);    
%        [LRImages W]=new_get_img_params(LRImages,5,128,128,64,64,model.SR);

    SR = imageToVector(model.SR);
            
    scgOptions = setupSCGOptions(solverParams);
    disp('Minimize objective function...');
    [SR, ~, flog, ~, ~]  = scg(@mapfunc, SR', scgOptions, @mapfunc_grad, model, LRImages, W, W');
    
    % Reshape parameter vector to a 2D image.
    SR = vectorToImage(SR, model.magFactor * lrDim);
    disp('DONE!');    
    
function scgOptions = setupSCGOptions(solverParams)

    scgOptions = zeros(1,18); 
    scgOptions(1) = 1;        
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
    r = LR - W*SR;
    dataTerm = sum((r.^2));
    
    % Evaluate image prior for regularization the super-resolved estimate.
    prior = model.imagePrior.function(SR, size(model.SR));
    
    % Calculate objective function.
    f = dataTerm + model.imagePrior.weight * prior;
                
function grad = mapfunc_grad(SR, model, LR, W, Wt)
    
    if ~iscolumn(SR)
        % Reshape to column vector. 
        SR = SR';
    end
    
    % Calculate gradient of the data fidelity term w.r.t. the
    % super-resolved image.
    r = LR - W*SR;
    dataTerm_grad = - 2 * (Wt * r);
    % Calculate gradient of the regularization term w.r.t. the 
    % super-resolved image.
    prior_grad = model.imagePrior.gradient(SR, size(model.SR));
   
    % Sum up to total gradient
    grad = dataTerm_grad + model.imagePrior.weight * prior_grad;
    grad = grad';