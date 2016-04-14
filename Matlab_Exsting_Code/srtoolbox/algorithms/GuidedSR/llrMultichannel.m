function [SR, report] = llrMultichannel(LRImages, model, varargin)
    
    % Get inter-channel prior parameters.
    if nargin > 2
        mu = varargin{1};
    else
        mu = 1.0;
    end
    if nargin > 3
        radius = varargin{2};
    else
        radius = 1;
    end
    if nargin > 4
        eps = varargin{3};
    else
        eps = 1e-4;
    end
    if nargin > 5
        isAdaptive = varargin{4};
    else
        isAdaptive = true;
    end
    
    if nargin > 6
        solverParams = varargin{5};
    else
        % Use default solver parameters
        solverParams = SRSolverParams;
    end
    
    % Get ground truth image for evaluation purposes if available.
    if nargin > 7
        gt = varargin{6};
    end
    
    % Assemble equation system to be solved consisting of low-resolution
    % observations and the system matrix.
    for c = 1:length(LRImages)  
        I_guidance = LRImages{c};
        lrDim{c} = size(I_guidance(:,:,1));
        if solverParams.verbose
            disp('Assemble equation system (system matrix)...');
        end
        [W{c}, LRImages{c}] = composeSREquationSystem(LRImages{c}, model);
        Wt{c} = W{c}';
    end
        
    if ~isfield(model, 'SR') || isempty(model.SR);
        % Compute initial guess for the super-resolved image. We use the
        % "average image" estimated from the system matrix and photometric
        % parameters (if available).
        if solverParams.verbose
            disp('Compute average image used as initial guess...');
        end
        for c = 1:length(LRImages)
            SR{c} = getInitialSRImage(W{c}, LRImages{c}, model.photometricParams);
        end
    else
        % Use initial guess provided by the user and reshape 2-D image into
        % parameter vector.
        for c = 1:length(model.SR)
            SR{c} = imageToVector(model.SR{c});
        end
    end
    
    scgOptions = setupSCGOptions(solverParams);
    if solverParams.verbose
        disp('Minimize objective function...');
    end
    
    for iter = 1:solverParams.maxIrlsIter
        
        if exist('gt', 'var')
            % Evaluate super-resolved data for the current iteration if
            % ground truth data is provided.
            for k = 1:length(LRImages)
                if iscell(gt)
                    I_gt(:,:,k) = gt{k};
                else
                    I_gt(:,:,k) = gt(:,:,k);
                end
                I_SR(:,:,k) = vectorToImage(SR{k}, model.magFactor * lrDim{k});
            end
            report.SR{iter} = I_SR;
            report.psnr_trace(iter) = psnr(I_SR, I_gt);
        end
        
        % -----------------------------------------------------------------
        % Step 1: Setup locally linear regression model
        % -----------------------------------------------------------------
        % Determine LLR filter coefficients and confidence weights for all
        % pairs of channels.
        for m = 1:length(LRImages)
            for n = 1:length(LRImages)
                if m ~= n
                    
                    % Determine the LLR filter coefficients for the current
                    % pair of channels
                    [A{m,n}, b{m,n}, R] = computeLLRFilterCoefficients(vectorToImage(SR{m}, model.magFactor * lrDim{m}), vectorToImage(SR{n}, model.magFactor * lrDim{m}), radius, eps);
                    % Determine adaptive confidence weights
                    if isAdaptive
                        % Select threshold to discriminate inliers and
                        % outliers in the LLR model.
                        thresh = 1.4826 * mad(R(:), 1);
                        % Compute confidence weights with selected
                        % threshold
                        w{m,n} = computeLLRConfidenceWeights( imageToVector(R), thresh );                                    
                    else
                        % Use uniform weights.
                        w{m,n} = ones(size(imageToVector(I_input)));
                    end
                    
                    % Reshape coefficients to matrix/vector notation.
                    A{m,n} = spdiags(imageToVector(A{m,n}), 0, numel(A{m,n}), numel(A{m,n}));
                    b{m,n} = imageToVector(b{m,n});
                end
            end
        end
                
        % -----------------------------------------------------------------
        % Step 2: Update multi-channel image
        % -----------------------------------------------------------------
        % Reshape cell arrays to parameter vector.
        SRvec = [];
        LRvec = [];
        for m = 1:length(SR)
            SRvec = [SRvec; SR{m}];
            LRvec = [LRvec; LRImages{m}];
        end
        % Perform SCG iterations to update the multi-channel image.
        SR_new = scg(@mapfunc, SRvec', scgOptions, @mapfunc_grad, model, LRvec, W, Wt, A, b, w, mu);
        
        % Reshape parameter vector into cell array.
        SR = mat2cell(SR_new', size(W{1}, 2) * ones(1,length(W)), 1);
               
    end
    
    % Reshape parameter vector to a 2D image.
    for c = 1:length(LRImages)
        SR{c} = vectorToImage(SR{c}, model.magFactor * lrDim{c});
    end
    if solverParams.verbose
        disp('DONE!');
    end
    
function scgOptions = setupSCGOptions(solverParams)

    scgOptions = zeros(1,18); 
    scgOptions(1) = solverParams.verbose;        
    scgOptions(2) = solverParams.tolX;   
    scgOptions(3) = solverParams.tolF;
    scgOptions(9) = solverParams.gradCheck;
    scgOptions(10) = solverParams.maxFunEvals;     
    scgOptions(14) = solverParams.maxIter;
    
function f = mapfunc(SR, model, LR, W, ~, A, b, w, mu)
    
    if ~iscolumn(SR)
        % Reshape to column vector. 
        SR = SR';
    end
    
    % Convert parameter vectors into cell array for further computations.
    SR = mat2cell(SR, size(W{1}, 2) * ones(1,length(W)), 1);
    LR = mat2cell(LR, size(W{1}, 1) * ones(1,length(W)), 1);
    
    % Compute the data fidelity term channel-wise.
    fData = 0;
    for m = 1:length(SR)
        model_m = model;
        if ~isempty(model.confidence)
            % Use confidence weights for the observations if these weights
            % are provided.
            model_m.confidence = model.confidence{m};
        else
            % Use uniform weights.
            model_m.confidence = ones(size(LR{m}));
        end
        fData = fData + mapDataTerm(SR{m}, model_m, LR{m}, W{m});
    end
    
    % Evaluate intra-channel prior for each channel.
    fIntra = 0;
    for m = 1:length(SR)
        fIntra = fIntra + model.imagePrior.weight(m) * model.imagePrior.function(SR{m}, model.imagePrior.parameters{1:end});
    end
    
    % Evaluate inter-channel prior for all pairs of channels.
    fInter = 0;
    for m = 1:length(W)
        for n = 1:length(W)
            if m ~= n
                fInter = fInter + mu(m,n) * sum( w{m,n} .* (A{m,n} * SR{m} + b{m,n} - SR{n}).^2 );
            end
        end
    end
        
    % Calculate objective function.
    f = fData + fIntra + fInter;
                
function grad = mapfunc_grad(SR, model, LR, W, Wt, A, b, w, mu)
    
    if ~iscolumn(SR)
        % Reshape to column vector. 
        SR = SR';
    end
    
    % Convert parameter vectors into cell array for further computations.
    SR = mat2cell(SR, size(W{1}, 2) * ones(1,length(W)), 1);
    LR = mat2cell(LR, size(W{1}, 1) * ones(1,length(W)), 1);
    
    % Calculate gradient of the data fidelity term w.r.t. the individual
    % channels.
    gData = [];
    for m = 1:length(W)
        model_m = model;
        if ~isempty(model.confidence)
            % Use confidence weights for the observations if these weights
            % are provided.
            model_m.confidence = model.confidence{m};
        else
            % Use uniform weights.
            model_m.confidence = ones(size(LR{m}));
        end
        gData = [gData; mapDataTerm_gradImage(SR{m}, model_m, LR{m}, W{m}, Wt{m})];
    end
    
    % Calculate gradient of the intra-channel prior.
    gIntra = [];
    for m = 1:length(W)
        gIntra = [gIntra; model.imagePrior.weight(m) * model.imagePrior.gradient(SR{m}, model.imagePrior.parameters{1:end})];
    end
    
    % Calculate gradient of the inter-channel prior for all pairs of
    % channels.
    gInter = [];
    for i = 1:length(W)
        g1 = 0;
        for n = 1:length(W)
            if n ~= i
                g1 = g1 + mu(i,n) * w{i,n} .* ( 2*A{i,n} * (A{i,n} * SR{i} + b{i,n} - SR{n}) );
            end
        end  
        g2 = 0;
        for m = 1:length(W)
            if m ~= i
                g2 = g2 - mu(m,i) * w{m,i} .* ( 2*(A{m,i} * SR{m} + b{m,i} - SR{i}) );
            end
        end      
        gInter = [gInter; g1 + g2];
    end
      
    % Compute overall gradient of the energy function.
    grad = gData + gIntra + gInter;
    grad = grad';
    