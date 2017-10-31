function [H, rand_colors, FigHandle] = visualizex(points, orders, ids, edges, name, logdir, showtext, showColor, nofigure)
    if nargin < 7
        showtext = true;
    end

    if nargin < 8
        showColor = false;
    end

    if nargin < 9
        nofigure  = false;
    end
    
    if numel(showColor) > 3
        rand_colors = showColor;
        showColor = true;
    else
        % random colors
        rand_colors = rand(500, 3);
    end

    % if not given orders
    if length(orders) == 0
        orders = 1: size(points, 2);
    end

    if ~nofigure
        FigHandle = figure;
        set(FigHandle, 'Position', [0, 0, 600, 1200]);
    end

    pointsm = [];
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
        if showtext
            if length(ids) > 0
                text(start(1), start(2), start(3), [num2str(ids(it))], 'FontSize', 14);
            end
        end
        pointsm = [pointsm, start];
    end

    if length(edges) > 0
        pointsm = [pointsm mean(pointsm, 2)];
        orders  = [orders; 0];
        for ed = 1: size(edges, 1)
            e = [pointsm(:, edges(ed, 1)), pointsm(:, edges(ed, 2))];
            plot3(e(1, :), e(2, :), e(3, :), '-', 'LineWidth', 2, 'Color', [1, 0, 0]); %axis equal; 
            hold on;
        end
        plot3(pointsm(1, end), pointsm(2, end), pointsm(3, end), 'Sk','MarkerFaceColor', 'k', 'MarkerSize', 5);
    end
    title(name); hold off;
    axis equal;
    %view(-6, -34)
    if length(logdir) > 0   
        savefig([logdir '/' name '.fig']);
    end
end