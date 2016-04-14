function prior = gaussianPrior_grad(x, imsize, h)

    X = vectorToImage(x, imsize);
    z = imageToVector( imfilter(X, h) );
    g = 2*z;
    
    G = vectorToImage(g, imsize);
    prior = imageToVector( imfilter(G, h) );