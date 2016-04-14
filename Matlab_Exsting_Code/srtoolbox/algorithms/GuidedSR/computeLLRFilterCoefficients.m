function [A, b, R] = computeLLRFilterCoefficients(I_guidance, I_input, radius, eps)
    
    % Determine the filter coefficients by guided filtering of the given
    % pair of channels.
    [A, b] = guidedfilter(I_guidance, I_input, radius, eps);
    
    % Get (filtered) residual error associated with the channels and 
    % filter coefficents.
    R =  A .* I_guidance + b - I_input;
    R = medfilt2(R, [3 3]);