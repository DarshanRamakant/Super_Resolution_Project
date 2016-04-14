function [f, g] = germanMcClureFunction(z, sigma)

    f = z.^2 ./ (z.^2 + sigma);
    if nargout > 1
        g = (2*sigma * z) ./ ((z.^2 + sigma).^2);
    end