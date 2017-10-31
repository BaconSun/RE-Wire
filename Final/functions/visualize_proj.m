function [H, rand_colors, FigHandle, I] = visualize_proj(points, orders, features, edges, ...
                                                         name, logdir, showtext, showColor, model)
    rand_colors = rand(500, 3);
    % if not given orders
    if length(orders) == 0
        orders = 1: size(points, 2);
    end

    FigHandle = figure;
    set(FigHandle, 'Position', [0, 0, 1000, 2000]);

    pointsm = [];
    subplot(221);
    for it  = 1: length(points)
        if showColor
            % scatter3(points{it}(1, :), points{it}(2, :), points{it}(3, :), 'Marker', '.', 'MarkerEdgeColor', rand_colors(orders(it) + 1, :), 'LineWidth', 2);  hold on;
            H(it) = plot3(points{it}(1, :), points{it}(2, :), points{it}(3, :), '-', 'Color', rand_colors(orders(it) + 1, :), 'LineWidth', 1.5, 'MarkerSize', 5);  hold on;
            set(H(it), 'UserData', orders(it));
        else
            H(it) = plot3(points{it}(1, :), points{it}(2, :), points{it}(3, :), '-', 'Color', rand_colors(1, :), 'LineWidth', 1.5, 'MarkerSize', 5);  hold on;
            set(H(it), 'UserData', orders(it));
        end

        start  = points{it}(:, round(size(points{it}, 2)/2)); % + [0.1, 0.1, 0.1]'
%         if length(features(it, :)) == 3
%             T(it) = text(start(1), start(2), start(3), [sprintf('%.3f', features(it, 1)) ',' sprintf('%.3f', features(it, 2)) ',' sprintf('%.3f', features(it, 3))], 'FontSize', 14);
%         else
%             T(it) = text(start(1), start(2), start(3), [sprintf('%.3f', features(it))], 'FontSize', 14);
%         end
%         if ~showtext
%             set(T(it), 'Visible', 'off');
%         end
        pointsm = [pointsm, start];
    end
    axis equal;

    % projections
    T1 = read_curvepoints([model 'curvepoints1.txt']);
    T2 = read_curvepoints([model 'curvepoints2.txt']);
    T3 = read_curvepoints([model 'curvepoints3.txt']);
    
    % camera parameters
    M1 = read_camera([model 'camera1.txt']);
    M2 = read_camera([model 'camera2.txt']);
    M3 = read_camera([model 'camera3.txt']);
    
    subplot(222); 
    plot(T1(:, 1),  -T1(:, 2),  '.r'); hold on;
    for it = 1: length(points)   
        V1 = project(M1, points{it});
        I1(it) = plot(V1(:, 1),  -V1(:, 2),  '.b'); hold on;
        set(I1(it), 'Visible', 'off');
        axis equal;
    end
    
    subplot(223);
    plot(T2(:, 1),  -T2(:, 2),  '.r'); hold on;
    for it = 1: length(points)   
        V2 = project(M2, points{it});
        I2(it) = plot(V2(:, 1),  -V2(:, 2),  '.b'); hold on;
        set(I2(it), 'Visible', 'off');
        axis equal;
    end
    
    subplot(224);
    plot(T3(:, 1),  -T3(:, 2),  '.r'); hold on;
    for it = 1: length(points)   
        V3 = project(M3, points{it});
        I3(it) = plot(V3(:, 1),  -V3(:, 2),  '.b'); hold on;
        set(I3(it), 'Visible', 'off');
        axis equal;
    end
    
    I{1} = I1;
    I{2} = I2;
    I{3} = I3;
    %title(name); hold off;
    %axis equal;
    %view(-6, -34)
    savefig([logdir '/' name '.fig']);
end

