function M = read_camera2(filename)
    matrix = importdata(filename);
    K  = matrix(4:6, 1:3);
    R  = matrix(7:9, 1:3);
    T  = matrix(10, 1:3)';

    RT = [R T];

    M  = K * RT;
end