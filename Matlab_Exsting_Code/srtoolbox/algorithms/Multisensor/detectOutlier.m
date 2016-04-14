function beta = detectOutlier(I_warped, I_fixed, refIdx, thresh, winSize, method, outSize)
    
    beta = [];
    for k = 1:size(I_warped, 3)  
        if k ~= refIdx          
            % Reorganize image pair to complex array for block processing.
            I = I_warped(:,:,k) + 1j*I_fixed(:,:,k);
            
            if strcmp(method, 'patch')
                % Evaluate patch-wise similarity. Resize confidence map to the
                % desired dimension if required.
                beta_k = imresize( blockproc(I, [winSize, winSize], @ncc), outSize );
            elseif strcmp(method, 'sliding')
                % Evaluate similarity in a sliding window scheme. Resize confidence map to the
                % desired dimension if required.
                beta_k = imresize( nlfilter(I, [winSize, winSize], @ncc), outSize );
            end
            
            % Reject confidence measures below the given threshold.
            beta_k((beta_k < thresh) | (isnan(beta_k))) = 0;
            % Truncate confidence map to 0...1.
            beta_k( beta_k > 1 ) = 1;
            
            % Linearize confidence map for further processing.
            beta_k = beta_k';
            beta = [beta; beta_k(:)]; %#ok<AGROW>          
        else
            % Set confidence map to the all-one vector for the reference 
            % frame.
            beta = [beta; ...
                    ones(outSize(1)*outSize(2), 1)]; %#ok<AGROW>
        end     
    end

function weight = ncc(block)
    
    if isstruct(block)
        block = block.data;
    end
    I1 = double( real(block) );
    I2 = double( imag(block) );
    weight = corr(I1(:), I2(:));
    weight = 0.5 * (1 + weight);