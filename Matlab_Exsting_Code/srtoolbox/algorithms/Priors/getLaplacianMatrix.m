%% GETLOGMATRIX - Determine matrix operator for Laplacian of Gaussian (LoG)
%
% D = GETLOGMATRIX(IMSIZE) determines the LoG matrix D for an image X of
% size IMSIZE such that z = D*X represents the LoG result of X.
function D = getLaplacianMatrix(imsize)

    % Call the MEX function to generate the LoG matrix.
    D = getLaplacianMatrix_mex(imsize(1), imsize(2));