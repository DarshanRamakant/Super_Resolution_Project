function f = huberPrior(x, imsize, h, thresh, weights)
    
    if nargin < 4
        % Default value for Huber function threshold parameter.
        thresh = 5e-3;
    end
    
    if nargin < 5
        % Use uniform weights for spatially adaptive version of the Huber
        % prior.
        weights = ones(size(x));
    end
    
    % Compute the high-pass filtered version of the image x using the
    % circular convolution matrix.
    z = imfilter(vectorToImage(x, imsize), h, 'symmetric');
    
    % Compute Huber loss function.
    huber = thresh * (sqrt(1 + (z./thresh).^2) - 1);
    %f = sum(weights .* huber);
    f = sum(weights .* imageToVector(huber));