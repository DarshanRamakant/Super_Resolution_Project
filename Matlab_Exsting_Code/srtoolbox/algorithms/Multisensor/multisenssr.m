function [SR, modelTransformed] = multisenssr(LRImages, model, varargin)
  
    if nargin > 2
        solverParams = varargin{end};
    else
        % Use default solver parameters
        solverParams = SRSolverParams;
    end
        
    % Transform motion from guidance images to low-resolution data.
    modelTransformed = model;
    numFrames = size(LRImages, 3);
    for k = 1:numFrames
        vx = model.motionParams{k}.vx;
        vy = model.motionParams{k}.vy;
        if any(vx) & any(vy) 
            [vx, vy] = projectMotion(vx, vy, [size(LRImages,1) size(LRImages,2)]);
        else
            refFrameIdx = k;
        end
        
        modelTransformed.motionParams{k}.vx = vx;
        modelTransformed.motionParams{k}.vy = vy;
    end
    
    if isempty(modelTransformed.photometricParams)
        % Derive photometric registration for low-resolution images from
        % transformed motion parameters.
        for k = 1:numFrames
            if k == refFrameIdx
                modelTransformed.photometricParams.mult(:,:,k) = 1;
                modelTransformed.photometricParams.add(:,:,k) = 0;
            else
                warp = modelTransformed.motionParams{k};
                [cMult, cAdd] = rangeCorr(medfilt2(LRImages(:,:,refFrameIdx), [3 3]), medfilt2(warpImage(LRImages(:,:,k), warp), [3 3]), 'msac', 0.025, 200);        
                modelTransformed.photometricParams.mult(:,:,k) = cMult;
                modelTransformed.photometricParams.add(:,:,k) = cAdd;
            end
        end
    end
    
    % Perform MAP super-resolution for low-resolution images with 
    % transformed motion.
    SR = superresolve(LRImages, modelTransformed, 'map', solverParams);