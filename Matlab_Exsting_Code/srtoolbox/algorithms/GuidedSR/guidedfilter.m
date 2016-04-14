function [A, b] = guidedfilter(I_guidance, I_input, radius, eps)
% O(1) time implementation of guided filter according to He et al.

    [h, w] = size(I_guidance);
    N = boxfilter(ones(h, w), radius);
    
    % Compute mean and covariances using the boxfilter.
    mean_guidance = boxfilter(I_guidance, radius) ./ N;
    mean_input = boxfilter(I_input, radius) ./ N;
    mean_Ip = boxfilter(I_guidance.*I_input, radius) ./ N;
    cov_Ip = mean_Ip - mean_guidance .* mean_input;
    mean_II = boxfilter(I_guidance.*I_guidance, radius) ./ N;
    var_I = mean_II - mean_guidance .* mean_guidance;

    % Compute coefficients for the patches.
    A = cov_Ip ./ (var_I + eps); 
    b = mean_input - A .* mean_guidance; 

    % Compute mean of all coefficients per patch using the boxfilter..
    A = boxfilter(A, radius) ./ N;
    b = boxfilter(b, radius) ./ N;