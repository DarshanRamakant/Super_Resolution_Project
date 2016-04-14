function g = lorentzianPrior_grad(x, imsize, h, alpha)
    
    % Compute filtered image.
    X = vectorToImage(x, imsize);
    z = imageToVector( imfilter(X, h) );
    
    % Compute the gradient of the Lorentzian function.
    lor_g = (2 * (z+eps)) ./ (2*alpha^2 + (z+eps).^2);
    
    % Compute gradient of the regularization term under the given high-pass
    % filter.
    lor_g = vectorToImage(lor_g, imsize);
    g = imageToVector( imfilter(lor_g, h) );
    
