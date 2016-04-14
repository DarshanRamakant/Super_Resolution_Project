function weights = getSpatiallyAdaptiveRegularizationWeights(I_input, I_guidance, patchSizeRadius, contrastFactor, thresh)
% Get spatially adaptive regularization weights for (adaptive) multi-sensor
% super-resolution.
    
    % Convert RGB to HSV space for spatially adaptive regularization.
    if ndims(I_guidance) == 3
        I_guidance = rgb2hsv(I_guidance);
        % Use V channel in the HSV color space.
        I_guidance = I_guidance(:,:,3);
    end
    
    % Edge detection on guidance image (i.e. photometric data)
    edgeMap = im2double(edge(I_guidance, 'sobel'));
    edgeMap(edgeMap > 1) = 1;
    edgeMap(edgeMap < 0) = 0;
    
    % Compute patch-wise similarity between input and guidance data (i.e. 
    % between range and photometric data) --> local mutual information.
    sim = computeLMI(I_guidance, imresize(I_input, size(I_guidance)), patchSizeRadius);
    % Thresholding for less significant correlation values (optional).
    if nargin > 4
        sim(sim < thresh) = 0;
    end
    
    % Assemble spatially adaptive regularization weights
    % Resample similarity to the dimension of the input data.
    sim = imresize(edgeMap .* sim, size(I_input));
    % Avoid negative similarity values.
    sim(sim < 0) = 0;
    % Compute spatially adaptive weights from local similarity.
    weights = exp(- sim ./ contrastFactor);
    weights = imageToVector(weights);