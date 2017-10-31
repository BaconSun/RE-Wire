function U = project(P, points)
    if size(points, 2) == 3
        if size(points, 1) ~= 3
            points = points';
        end
    end
    if size(points, 1) == 3
        points = [points; ones(1, size(points, 2))];
    else if size(points, 1) == 7
        points = points(1:3, :);
        points = [points; ones(1, size(points, 2))];
    
        end
    end
    if isempty(points)
        U =[];
    else
        U  = P * points;
        U  = U ./ repmat(U(3, :), 3, 1);
        U  = round(U(1:2, :))';   % get image
    end
end
