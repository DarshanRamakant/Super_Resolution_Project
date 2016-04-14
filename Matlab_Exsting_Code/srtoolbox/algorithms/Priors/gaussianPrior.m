function prior = gaussianPrior(x, imsize, h)
    
    % High-pass filtering of input image.
    X = vectorToImage(x, imsize);
    z = imageToVector( imfilter(X, h) );
    
    % Apply L2 loss function.
    prior = sum(z.^2);