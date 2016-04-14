function imSeq = readRangeSequence(dirName, pattern, roi)

    % Get image filenames from directory.
    filenames = dir( [dirName, '\', pattern] );
    numFrames = length(filenames);
    
    % Read first image to allocate space for image sequence. All images
    % must be of the same dimension.
    firstImg = dlmread( [dirName, '\', filenames(1).name] );
    imSeq = zeros( size(firstImg, 1), size(firstImg, 2), numFrames );
    imSeq(:,:,1) = firstImg;
    
    % Read other images of the sequence.
    for k = 2:numFrames
        imSeq(:,:,k) = dlmread( [dirName, '\', filenames(k).name] );
    end
    
    if nargin > 2 && (~isempty(roi))
        % Crop images according to specified ROI.
        rowStart = roi(1);
        colStart = roi(2);
        h = roi(3);
        w = roi(4);
        imSeq = imSeq(rowStart:(rowStart + h - 1), colStart:(colStart + w - 1), :);
    end