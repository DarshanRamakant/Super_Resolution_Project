function lmi = computeLMI(I_guidance, I_input, patchSizeRadius)
% Compute local mutual information with given patch size as similarity
% measure between an guidance image (RGB image) and an input image (range
% image).

    dimY = size(I_input, 1);
    dimX = size(I_input, 2);
    lmi = zeros(size(I_input));
    for y = 1:size(I_input, 1)
        for x = 1:size(I_input, 2);
            
            % Extract neighbourhood for current position.
            y0 = max((y-patchSizeRadius), 1);
            y1 = min((y+patchSizeRadius), dimY);
            x0 = max((x-patchSizeRadius), 1);
            x1 = min((x+patchSizeRadius), dimX);
            I = I_guidance(y0:y1, x0:x1);
            J = I_input(y0:y1, x0:x1);
            
            % Compute mutual information and joint entropy for this window.
            % We use the MutualInfo package written by Hanchuan Peng to
            % compute both quantities.
            mi = mutualinfo(uint8(255*I(:)), uint8(255*J(:)));
            h = jointentropy(uint8(255*I(:)), uint8(255*J(:)));
            
            % Compute local mutual information (LMI).
            lmi(y, x) = mi / h; 
        end
    end