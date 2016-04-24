%function prior = gaussianPrior(x, imsize, h)
function prior = gaussianPrior(x, imsize)
    whos x;
    % High-pass filtering of input image.
    X = vectorToImage(x, imsize);
   
    %h = fspecial('gaussian', 3, 0.5);
    h = [1 1 1;
     1 -8 1;
     1  1 1];
    z = imageToVector( imfilter(X, h) );
      
    %z = imageToVector( imgaussfilt(X,0.5));

    
    % Apply L2 loss function.
    prior = sum(z.^2);