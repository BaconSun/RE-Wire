%unction cleanTSP(task, verbose)
    verbose=true
    task = setup('turtle_real')
    colors  = [255 0 0; 0 0 255; 0 255 255; 255 255 0; 0 255 0; 255 0 255]/255;

    %task   = setup('cart');
    model  = task.model;
    modelpath = ['./tasks/' model '/'];
    mu     = task.select(1);
    lambda = task.select(2);
    tsplam = task.select(3);
    suffix = task.suffix;

    % check if selected
    if ~exist([modelpath 'slog' suffix '/mu=' num2str(mu) '_l=' num2str(lambda)])
        error(['there is no selected curves for mu=', num2str(mu)]);
    end

    disp(['start TSP: mu=' num2str(mu) ', lambda=' num2str(lambda)]);
    load([modelpath 'slog' suffix '/mu=' num2str(mu) '_l=' num2str(lambda) '/data.mat']);
    if task.type == 'm'
        cornerfile = [modelpath 'slog/mu=' num2str(mu) '_l=' num2str(lambda) '/corners.mat'];
        if exist(cornerfile)
            load(cornerfile);
        else
            M = [];
        end
    end

    % TSP logs
    if ~exist([modelpath 'log'])
        mkdir([modelpath 'log']);
    end

    logdir=[modelpath 'log/mu=' num2str(mu) '_l=' num2str(lambda)];
    if ~exist(logdir)
        mkdir(logdir);
    end

    load([logdir '/tsp.mat'])
    disp('------------------------------------------------------------------')
    
    points0 = points;
    tag = ones(L-1, 1);

    start_points = find(orders==1)
    end_point = max(max(edges))

    for num = 1: length(start_points)
        a0 = start_points(num);
        pathx = [a0];
        a1  = a0;

        while true
            s = find(edges(:, 1)==a1);
            if edges(s, 2) == end_point
                break
            end
            a1 = edges(s, 2);
            pathx = [pathx a1];
        end
        path_all{num} = pathx;
    end
    
    for num = 1: length(path_all)

        path1 = path_all{num};
        path2 = [path1(1: end-1); path1(2: end)]';
        L     = numel(path1);

        direction    = zeros(size(path1));
        for it = 1: length(path1) - 1
            a = path1(it);
            b = path1(it + 1);
            d1 = D1(a, b); d2 = D2(a, b); d3 = D3(a, b); d4 = D4(a, b);
            if it == 1
                direction(it) = (min(d4, d1) > min(d2, d3));
            end
            direction(it + 1) = (min(d1, d3) > min(d2, d4));
        end


        % filter-1
        a_tags  = ones(L, 1);
        a_tags(1) = 0;
        a_tags(end) = 0;
        a_tags(end-1) = 0;

        a1_tags = a_tags;
        a2_tags = a_tags;

        % fiter-2
        points = {};
        for k = 1: L
            points{k} = points0{path1(k)};
        end
        

        for k =1 : (L-3)
            A = points{k}; B = points{k+1}; C = points{k+2}; D = points{k+3}; 

            if direction(k) == 0
                Aend = points{k}(:, 1);
            else
                Aend = points{k}(:, end);
            end

            if direction(k+1) == 0
                Bhead = points{k+1}(:, end);
                Bend  = points{k+1}(:, 1);
            else
                Bhead = points{k+1}(:, 1);
                Bend  = points{k+1}(:, end);
            end

            if direction(k+2) == 0
                Chead = points{k+2}(:, end);
                Cend  = points{k+2}(:, 1);
            else
                Chead = points{k+2}(:, 1);
                Cend  = points{k+2}(:, end);
            end
            k
            if direction(k+3) == 0
                Dhead = points{k+3}(:, end);
            else
                Dhead = points{k+3}(:, 1);
            end

            % condition-1
            dda = norm(Dhead - Aend)      
            dcb = norm(Chead - Bend)
            if (dda < 0.08) & (dcb < 0.02)
                a_tags(k+1)  = a_tags(k+1) * 1;
                a1_tags(k+1) = a1_tags(k+1) * 1;
            else
                a_tags(k+1)  = a_tags(k+1) * 0;
                a1_tags(k+1) = a1_tags(k+1) * 0;
            end

            % condition-2
            V1 = normalize(Bend - Bhead);
            V2 = normalize(Chead - Cend);
            av1v2 = dot(V1, V2)
            if av1v2 > 0.9
                a_tags(k+1)  = a_tags(k+1) * 1;
                a2_tags(k+1) = a2_tags(k+1) * 1;
            else
                a_tags(k+1)  = a_tags(k+1) * 0;
                a2_tags(k+1) = a2_tags(k+1) * 0;
            end



            % condition-3
            dbb = norm(Bend - Bhead)
            dcc = norm(Cend - Chead)

        end

        if verbose
            visualizex(points, a1_tags, 1:length(points), [], 'condition-1', logdir, true, colors);
            visualizex(points, a2_tags, 1:length(points), [], 'condition-2', logdir, true, colors);
            visualizex(points, a_tags,  1:length(points), [], 'TSP edges', logdir, true, colors);
        end
    end

    return

    a_tags = (a_tags + [0; a_tags(1:end-1)] > 0);

    pointss = [];
    pointse = [];
    points_ = {};
    for sk = 1:length(a_tags)
        if a_tags(sk) == 0
            points_{end+1} = points{sk};
            pointss = [pointss, points{sk}(:, 1)];
            pointse = [pointse, points{sk}(:, end)];
        end
    end

    %% use curve distance
    D1 = pdist2(pointss', pointse');
    D2 = pdist2(pointse', pointss');
    D3 = pdist2(pointse', pointse');
    D4 = pdist2(pointss', pointss');


    points = points_;
    LLL    = length(points);
    orders = 1:LLL;
    orders = orders';
    edges  = [LLL+1 1; orders(1:end-1) orders(2:end); LLL LLL+1];
    
    if verbose
        visualizex(points, orders, orders, edges, 'new TSP edges', logdir, true, true);
    end
    %visualize2(points, 1:length(points), [], 'TSP edges', logdir, true, true);
    %visualize2(points, d_tags, [], 'TSP edges', logdir, true, true);
    %visualize2(points, 1 - (1 - a_tags) .* (1 - d_tags), [], 'TSP edges', logdir, true, true);

    disp('save TSP infomation')

    % save([logdir '/tsp.mat'], 'edges', 'orders', 'points', 'D1', 'D2', 'D3', 'D4')
%end