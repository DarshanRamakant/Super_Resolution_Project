function [vx, vy] = projectMotion(vxGuidance, vyGuidance, outSize)
    
    if outSize == size(vxGuidance(:,:,1))
        vx = vxGuidance;
        vy = vyGuidance;
        return;
    end
    
    % Define function to transform displacement vector fields: Median
    % displacement in a local patch.
    localMedian = @(block_struct) ... 
                    median(double(block_struct.data(:)));
    
    % Determine patch size and associated scaling.
    patchSize = round(size(vxGuidance) ./ outSize);
    scale = 1 ./ patchSize;
    
    % Transform and rescale displacements in x-direction.
    vx = blockproc(vxGuidance, patchSize, localMedian);
    vx = scale(2) * vx;
    
    % Transform and rescale displacements in y-direction.
    vy = blockproc(vyGuidance, patchSize, localMedian);
    vy = scale(1) * vy;
    
    % Resize displacement vector fields if the size does not match with the
    % desired size (in case of non-integer patch size).
    if (size(vx,1) ~= outSize(1)) || (size(vx,2) ~= outSize(2)) 
        vx = imresize(vx, outSize);
        vy = imresize(vy, outSize);
    end