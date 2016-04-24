function SR = getInitialSRImage(W, LRImages, photoParams)
% GETINITIALSRIMAGE Initial guess for super-resolved image
%   GETINITIALSRIMAGE computes a rough estimate for the super-resolved
%   image from given model parameters.
%
%   X = GETINITIALSRIMAGE(W, Y, PHOTOPARAMS) computes the "average image" X
%   as initial guess based on the system matrix W, the low-resolution data 
%   Y and the photometric parameters PHOTOPARAMS.
%
%   See David P. Capel, Image Mosaicing and Super-resolution, PhD thesis,
%   2001 for details.

        % Compute "average image" without photometric parameters
        SR = (W' * LRImages) ./ (sum(W, 1)' + 1e-4);
    end