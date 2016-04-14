% GENERATESYNTHETICIMAGESEQUENCE - Generate synthetic image sequence.
%
% GENERATESYNTHETICIMAGESEQUENCE(GT, IMAGPARAMS, MOTIONLB, MOTIONUP,
% NOISESTD, NUMFRAMES) generates NUMFRAMES synthetic images I for ground
% truth image GT with uniform distributed random motion out of the range
% [MOTIONLB, MOTIONUP]. The random motion is returned in MOTIONPARAMS.
%
% See also: generateSyntheticImage
function [I, motionParams] = generateSyntheticImageSequence(gt, magFactor, psfWidth, motionLowerBound, motionUpperBound, noiseStd, numFrames)

%I = zeros(imagingParams.lrDim(1), imagingParams.lrDim(2), numFrames);
for k = 1:numFrames
    
    % Generate some (uniform distributed) random motion parameters.
    motionParams{k} = zeros(3, 3);
    if k > 1
        % Select uniform distributed motion parameters.
        rotAngle = motionLowerBound.rotAngle + (motionUpperBound.rotAngle - motionLowerBound.rotAngle) * rand(1);
        translation = motionLowerBound.translation + (motionUpperBound.translation - motionLowerBound.translation) .* rand(2, 1);
        A = [cosd(rotAngle) -sind(rotAngle) translation(1); ...
            sind(rotAngle) cosd(rotAngle)  translation(2); ...
            0              0               1 ];
        motionParams{k} = A;
    else
        % There should be no motion between the first image and the
        % ground truth.
        motionParams{k} = [1 0 0; ...
            0 1 0; ...
            0 0 1];
    end
    
    % Generate a new frame for the sequence with the given motion parameters.
    if size(gt, 3) == 1
        I(:,:,k) = generateSyntheticImage(gt, magFactor, psfWidth, motionParams{k}, noiseStd);
    else
        I(:,:,:,k) = generateSyntheticImage(gt, magFactor, psfWidth, motionParams{k}, noiseStd);
    end
end

