function visualize(points, orders, edges, name, logdir, showtext)
    if nargin < 6
        showtext = true;
    end

    % if not given orders
    if length(orders) == 0
        orders = 1: size(points, 2);
    end

    FigHandle = figure;
    set(FigHandle, 'Position', [100, 100, 800, 600]);

    pointsm = [];
    for it  = 1: length(points)
        plot3(points{it}(1, :), points{it}(2, :), points{it}(3, :));  hold on;
        start  = points{it}(:, round(size(points{it}, 2)/2)); % + [0.1, 0.1, 0.1]'
        if showtext
            text(start(1), start(2), start(3), [num2str(orders(it))], 'FontSize', 14);
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
    savefig([logdir '/' name '.fig']);
end