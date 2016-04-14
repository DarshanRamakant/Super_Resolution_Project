function f = lorentzianPrior(x, imsize, h, alpha)
    
    % Compute filtered image.
    X = vectorToImage(x, imsize);
    z = imageToVector( imfilter(X, h) );
    
    % Compute function value of the Lorentzian function.
    z = z + eps;
    lor = log(1 + 0.5*((z+eps) ./ alpha).^2);
    
    % Compute the objective value of the regularization term.
    f = sum(lor);
