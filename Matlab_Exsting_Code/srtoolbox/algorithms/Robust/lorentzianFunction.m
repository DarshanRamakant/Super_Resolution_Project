function [f, g] = lorentzianFunction(z, sigma)
    
    z = z + eps;
    f = log(1 + 0.5*(z ./ sigma).^2);
    if nargout > 1
        g = (2 * z) ./ (2*sigma^2 + z.^2);
    end