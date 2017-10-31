function pipe = merge_nodes(pipe_in)

    disp('merge nodes');

    % pipe in -----------------------
    clusters   = pipe_in.clusters;
    check      = pipe_in.check;
    points0    = pipe_in.points;
    curveNums0 = pipe_in.curveNums;
    model      = pipe_in.model;

    disp(['before: ' num2str(length(points0))]);


    %% merge clusters --------------------------------------------
    % check single nodes in each clusters
    temp     = [diff(clusters); 1];
    merges   = [temp(1); temp(1: end-1) .* temp(2: end)];

    % check single nodes in the same combined set.
    merges   = merges .* check;
    new_node = ([merges(1); diff(merges)] > 0) + (merges == 0);

    %[clusters check new_node merges]


    % get new clusters and new check
    new_clusters = clusters((clusters .* new_node) > 0);
    new_clusters = cumsum([new_clusters(1); diff(new_clusters)] > 0);
    new_nums = [1];
    for it = 2: length(new_clusters)
        if new_clusters(it) == new_clusters(it-1)
            new_nums(end) = new_nums(end) + 1;
        else
            new_nums(end+1) = 1;
        end
    end
    new_check    = check((check .* new_node) > 0);
    new_combined = [1; diff(new_check)] .* (1: length(new_check))';
    new_combined = new_combined(new_combined>0)-1;

    %[new_clusters new_check]
    %return

    % save new .off files:
    if ~exist([model 'combined_off'])
        mkdir([model 'combined_off']);
    end
    source_path  = [model 'combined_off/'];

    new_id       = 1;
    temp_m       = 0;
    temp_pts     = points0{1};
    new_F{1}     = ['new_node_' num2str(new_id) '.off'];

    for it = 2: length(new_node)
        if new_node(it) > 0
            write_off([source_path new_F{end}], temp_pts, int16.empty(0,3));    
            new_id       = new_id + 1; 
            new_F{end+1} = ['new_node_' num2str(new_id) '.off'];
            temp_pts     = [];
        end
        temp_pts = [temp_pts points0{it}];
    end

    if length(temp_pts) > 0
        write_off([source_path new_F{end}], temp_pts, int16.empty(0,3));  
    end


    % get visualize
    endDireNum = 5;

    pointss = [];
    pointse = [];
    for it = 1: length(new_F)   
        points1{it} = read_off([source_path new_F{it}]);
        pointss = [pointss, points1{it}(:, 1)];
        pointse = [pointse, points1{it}(:, end)];
        curveNums1(it) = size(points1{it}, 2);
    end
    endPoints(:, 1, :) = pointss;
    endPoints(:, 2, :) = pointse;

    %if pipe_in.visualize
    %    logdir = pipe_in.logdir;

        % visualize2(points0, clusters, [], 'Clusters', logdir, true, false);
        % visualize2(points0, check, [],    'Combined', logdir, false, true);
        % visualize2(points1, new_clusters, [], 'after merge: Clusters', logdir, false, true);
        % visualize2(points1, new_check, [],    'after merge: Combined', logdir, false, true);
    %end

    %% compute new C, E, D --------------------------------------------------

    % compute D & E
    c1 = [1 1 2 2]'; c2 = [1 2 1 2]';
    for x = 1: length(new_F)
        for y = 1: length(new_F)
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
            
            dotProduct      = (abs(gapDire' * endDire1) + abs(gapDire' * endDire2)) / 2;
            E(x, y)         = exp(-dotProduct);
        end
    end

    % compute new C (use original "C")
    oldC = (importdata([model 'outC.txt']) .* curveNums0' * 3).^2;
    C    = []; 

    for it = 1: length(new_node)
        if new_node(it) > 0
            C(end+1) = oldC(it);
        else
            C(end) = C(end) + oldC(it);
        end
    end
    C    = sqrt(C) ./ (curveNums1 * 3);

    % output pipe
    pipe.clusters = new_clusters;
    pipe.combined = new_combined;
    pipe.nums     = new_nums;
    pipe.path     = source_path;
    pipe.C        = C;
    pipe.D        = D;
    pipe.E        = E;
    pipe.F        = new_F;

    disp(['after: ' num2str(length(points1))]);
end




