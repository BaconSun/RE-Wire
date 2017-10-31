function V = read_curvepoints(filename)
    t1 = importdata(filename);
    t1 = reshape(t1(2: end), 5, []);
    t1 = t1(2:3, :);
    V  = [t1(2, :); t1(1, :)]';
end