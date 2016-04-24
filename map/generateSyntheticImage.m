% GENERATESYNTHETICIMAGE - Generate shifted, blurred downsampled and noisy
%                          images.
%
% I = GENERATESYNTHETICIMAGE(GT, IMAGINGPARAMS, MOTIONPARAMS, NOISESTD)
% gernates a shifted, downsampled blurred and noisy version of a given
% ground truth image GT.
% The imaging parameters (PSF and downsampling factors) and the noise level
% of an additive white Gaussian noise are specified by IMAGINGPARAMS and
% NOISESTD respectively. The motion relative to the ground truth image is
% specified by the motion parameters MOTIONPARAMS.
%
% See also: generateSyntheticImageSequence
function [I] = generateSyntheticImage(gt, magFactor, psfWidth, motionParams, noiseStd)

% Convert to double for further processing.
gt = im2double(gt);

% Get system matrix for given imaging and motion parameters.
if isempty(motionParams)
    motionParams = eye(3);
end
W = composeSystemMatrix(size(gt(:,:,1)) ./ magFactor, magFactor, psfWidth, motionParams);

for c = 1:size(gt, 3)

    % Transform 2D image to a 1D data vector.
    x = imageToVector(gt(:,:,c));
    y = W * x;  % Determine shifted, blurred and downsampled data.

    % Add Gaussian noise with specified standard deviation.
    y = y + noiseStd * randn(size(y));
	
	% Clipping
	y(y < 0) = 0;
	y(y > 1) = 1;

    % Reshape data vector to 2D image.
    I(:,:,c) = vectorToImage(y, size(gt(:,:,1)) ./ magFactor);
    
end