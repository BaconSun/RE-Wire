function vis_view(V, U, fig, r, g, b)
    BG     = V;
    U      = U' ;
    if length(size(BG)) == 2
        BG = repmat(BG, 1, 1, 3);
    end

    for it = 1: size(U, 2)
        BG(U(2, it) - 5: U(2, it) + 5, U(1, it) - 5: U(1, it) + 5, 1) = r;
        BG(U(2, it) - 5: U(2, it) + 5, U(1, it) - 5: U(1, it) + 5, 2) = g;
        BG(U(2, it) - 5: U(2, it) + 5, U(1, it) - 5: U(1, it) + 5, 3) = b;
    end
    figure;
    % imshow(BG)
    imwrite(BG, fig)
end
