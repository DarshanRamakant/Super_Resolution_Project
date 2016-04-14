%% Setup path
pathOfThisFile = mfilename('fullpath');
rootDir = fileparts(pathOfThisFile);
if ~isempty(rootDir)
    addpath(genpath(rootDir));
    disp('Super-Resolution Toolbox path setup has been successful.');
end

%% Run MEX compiler
disp('Compile MEX files...');
disp('...getLaplacianMatrix_mex');

currentPath = pwd;
cd([rootDir, '/algorithms/Priors']);
mex getLaplacianMatrix_mex.cpp -largeArrayDims;

disp('...composeSystemMatrix_mex');
cd([rootDir, '/algorithms/MAP']);
mex composeSystemMatrix_mex.cpp -largeArrayDims;
cd(currentPath);

disp('DONE!');