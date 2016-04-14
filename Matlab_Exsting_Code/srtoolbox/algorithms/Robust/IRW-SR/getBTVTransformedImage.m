function btvTransformedImage = getBTVTransformedImage(X, P, alpha0)

    % Pad image at the border to perform shift operations.
    Xpad = padarray(X, [P P], 'symmetric');
    % Consider shifts in the interval [-P, +P].
    btvTransformedImage = [];
    for l = -P:P
        for m = -P:P
            if l ~= 0 || m ~= 0
                % Shift by l and m pixels.
                Xshift = Xpad((1+P-l):(end-P-l), (1+P-m):(end-P-m));
                btvTransformedImage = [ btvTransformedImage; ...
                                        imageToVector( medfilt2( alpha0^(abs(l) + abs(m)) * (Xshift - X), [3, 3]) ) ];
            end
        end
    end

end