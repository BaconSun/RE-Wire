function [E, D] = recompute(source_path, new_F)
    endDireNum = 4;

    pointss = [];
    pointse = [];
    for it = 1: length(new_F)   
        points1{it} = read_coff([source_path new_F{it}]);
        pointss = [pointss, points1{it}(:, 1)];
        pointse = [pointse, points1{it}(:, end)];
        curveNums1(it) = size(points1{it}, 2);
    end
    endPoints(:, 1, :) = pointss;
    endPoints(:, 2, :) = pointse;


    %% compute new C, E, D --------------------------------------------------

    % compute D & E
    c1 = [1 1 2 2]'; c2 = [1 2 1 2]';
    for x = 1: length(new_F)
        for y = 1: length(new_F)
            
            if x == y
                E(x, y) = 0;
                D(x, y) = 0;
                continue;
            end

            Gap(:, :, x, y) = endPoints(:, c1, x) - endPoints(:, c2, y);
            
            D0(:, x, y)     = sqrt(sum(Gap(:, :, x, y) .^ 2, 1));
            [Dist, Idx]     = min(D0(:, x, y));
            D(x, y)         = Dist;       % -------- D
            
            gapDire         = normalize(Gap(:, Idx, x, y));


            if c1(Idx) == 1
                endDire1 = normalize(points1{x}(:, 5) - points1{x}(:, 1));
            else
                endDire1 = normalize(points1{x}(:, end-4) - points1{x}(:, end));
            end
            if c2(Idx) == 1
                endDire2 = normalize(points1{y}(:, 5) - points1{y}(:, 1));
            else
                endDire2 = normalize(points1{y}(:, end-4) - points1{y}(:, end));
            end

            %% ATTENTION !!!!
            %% There is a hidden BUG when computing E !!!            
            if sum(gapDire) == 0        
                gapDire = -endDire2 + endDire1;
                if sum(gapDire) == 0
                    error('ZERO')
                end
                gapDire = normalize(gapDire);
            end
            
            dotProduct  = (abs(gapDire' * endDire1) + abs(gapDire' * endDire2)) / 2;
            
            E(x, y)         = exp(-dotProduct);
        end
    end
end