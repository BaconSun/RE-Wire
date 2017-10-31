function N=normal1(Vertices)

    x  = Vertices(:, 1);
    y  = Vertices(:, 2);
    dy = gradient(y);
    dx = gradient(x);
    dm = sqrt(dx .^ 2 + dy .^ 2);
    dx = dx ./ dm;
    dy = dy ./ dm;

    N  = [-dy dx];
end