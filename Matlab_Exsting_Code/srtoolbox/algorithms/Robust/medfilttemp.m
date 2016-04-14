function SR_med = medfilttemp(LRImages, motionParams)
        
    for k = 1:size(LRImages,3)
        H = motionParams{k};
        Iwarped(:,:,k) = imtransform(LRImages(:,:,k), maketform('projective', H'), 'XData', [1 size(LRImages,2)], 'YData', [1 size(LRImages,1)]);
    end
    
    SR_med = zeros(size(LRImages,1), size(LRImages,2));
    for y = 1:size(SR_med,1)
        for x = 1:size(SR_med,2)
            I = Iwarped(y, x, :);
            SR_med(y,x) = median(I(I > 0));
        end
    end
    SR_med(isnan(SR_med)) = 0;