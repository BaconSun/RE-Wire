function M = read_camera(filename)
    matrix = importdata(filename);
    F  = matrix(1);
    cx = matrix(2);
    cy = matrix(3);
    T  = matrix(4:6);
    R  = reshape(matrix(7:end), 3, 3)';

    K  = [F 0  cx; 0 F cy; 0 0 1];
    RT = [R T];

    M  = K * RT;
end