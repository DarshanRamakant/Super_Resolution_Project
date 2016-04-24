%function prior = gaussianPrior_grad(x, imsize, h)
function prior = gaussianPrior_grad(x, imsize)

    X = vectorToImage(x, imsize);

    %h = fspecial('gaussian', 3, 0.5);
h = [1 1 1;
     1 -8 1;
     1  1 1];
    z = imageToVector( imfilter(X, h) );
    %z = imageToVector( imgaussfilt(X,0.5));

    g = 2*z;
    
    G = vectorToImage(g, imsize);
    prior = imageToVector( imfilter(G, h) );
    %prior = imageToVector( imgaussfilt(G,0.5));
