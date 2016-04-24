addpath('..');
addpath('../../synthetic');
addpath('../../algorithms');
addpath('../../algorithms/Priors');
addpath('../../algorithms/MAP');
addpath('../../common');
addpath('../../common/3rdparty/netlab3_3');

%% Get ground truth image
%gt = im2double( imread('cameraman.tif') );
gt = im2double( imread('lena.jpg') );

%% Generate low-resolution frames
magFactor = 2;      % Magnification factor
psfWidth = 0.4;     % Width of isotropic Gaussian PSF
noiseStd = 0.025;   % Standard deviation of additive Gaussian noise
numFrames = 16;     % Number of low-resolution frames
% Uniformly random distributed translation and rotation
motionLowerBound.translation = [-2; -2];    % Lower bound for translation
motionLowerBound.rotAngle = -0.5;           % Lower bound for rotation angle
motionUpperBound.translation = [2; 2];      % Upper bound for translation
motionUpperBound.rotAngle = 0.5;            % Upper bound for rotation angle
[LRImages, gtMotion] = generateSyntheticImageSequence(gt, magFactor, psfWidth, motionLowerBound, motionUpperBound, noiseStd, numFrames);

%% Super-resolution reconstruction
% Setup model parameters
model = SRModel;

% Set the desired magnification factor.
model.magFactor = magFactor;
% Set motion parameters used for super-resolution.
model.motionParams = gtMotion;

% Setup image prior using bilateral total variation (BTV) model
lambda = 0.0025;    % Regularization weight
P = 2;              % Window size
alpha = 0.6;        % Scale factor
%model.imagePrior = SRPrior('function', @btvPrior, 'gradient', @btvPrior_grad, 'weight', lambda, 'parameters', {magFactor * size(LRImages(:,:,1)), P, alpha});
% model.imagePrior = SRPrior('function', @gaussianPrior, 'gradient', @gaussianPrior_grad, 'weight', lambda, 'parameters', {magFactor * size(LRImages(:,:,1))});
srprior = struct('function',   @huberPrior, ...        % Function handle to prior (default: Huber)
              'gradient',   @huberPrior_grad, ...   % Function handle to gradient
              'weight',     0, ...                  % Regularization weight
              'parameters', [], ...                 % Additional parameters passed to prior function/gradient
              'weightEvaluationFunction', []);    
 SRprior = struct('function',   @gaussianPrior, ...        % Function handle to prior (default: Huber)
                     'gradient',   @gaussianPrior_grad, ...   % Function handle to gradient
                     'weight',     lambda, ...                  % Regularization weight
                     'parameters', {magFactor * size(LRImages(:,:,1))});    % Additional parameters passed to prior function/gradient
          
          %model.imagePrior = SRPrior('function', @gaussianPrior, 'gradient', @gaussianPrior_grad, 'weight', lambda, 'parameters', {magFactor * size(LRImages(:,:,1))});

          model.imagePrior = SRprior;
model.SR = imresize(LRImages(:,:,1), model.magFactor);

solverParams = SRSolverParams('maxFunEvals', 25, 'maxIter', 25, 'tolX', 1e-2, 'tolF', 1e-2, 'verbose', true);
%I_sr = superresolve(LRImages, model, 'map', solverParams);
 I_sr = mapsr(LRImages, model,solverParams);


figure; imshow([model.SR I_sr gt]);