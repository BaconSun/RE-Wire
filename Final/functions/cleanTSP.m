function cleanTSP(task, verbose)
    colors  = [255 0 0; 0 0 255; 0 255 255; 255 255 0; 0 255 0; 255 0 255]/255;

    %task   = setup('cart');
    model  = task.model;
    modelpath = ['./tasks/' model '/'];
    lambda = task.lambda(1);
    mu     = task.lambda(2);
    tsplam = task.lambda(3);
    suffix = task.suffix;

    % check if selected
    if ~exist([modelpath 'slog' suffix '/lambda=' num2str(lambda) '_mu=' num2str(mu)])
        error(['there is no selected curves for lambda=', num2str(lambda)]);
    end
    load([modelpath 'slog' suffix '/lambda=' num2str(lambda) '_mu=' num2str(mu) '/data.mat']);
    if task.type == 'm'
        cornerfile = [modelpath 'slog/lambda=' num2str(lambda) '_mu=' num2str(mu) '/corners.mat'];
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

    logdir=[modelpath 'log/lambda=' num2str(lambda) '_mu=' num2str(mu)];
    if ~exist(logdir)
        mkdir(logdir);
    end

    load([logdir '/tsp.mat'])
    disp('------------------------------------------------------------------')
    
    L   = length(points);
    tag = ones(L-1, 1);
    path1 = arrayfun(@(x)find(orders==x), 1: L);
    path2 = [path1(1: end-1); path1(2: end)]';

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
    points_ = {};
    for k = 1: L
        points_{orders(k)} = points{k};
    end
    points = points_;

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
   

    disp('save TSP infomation')

end