function prior = bepPrior(x, imsize, P, alpha, bepThreshold)

    if nargin < 3
        P = 1;
    end
    if nargin < 4
        alpha = 0.7;
    end
    
    % Reshape SR vector to image for further processing.
    X = vectorToImage(x, imsize);
    
    % Pad image at the border to perform shift operations.
    Xpad = padarray(X, [P P], 'symmetric');

    % Consider shifts in the interval [-P, +P].
    prior = 0;
    for l=-P:P
        for m=-P:P
            % Shift by l and m pixels.
            Xshift = Xpad((1+P-l):(end-P-l), (1+P-m):(end-P-m));
            e = Xshift(:) - X(:);
            prior = prior + alpha.^(abs(l)+abs(m)) .* sum( bepThreshold * sqrt(e.^2 + bepThreshold^2) );      
        end
    end
    
    
    

