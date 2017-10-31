function runTSP(task, given_m, verbose)

    model  = task.model;
    modelpath = ['./tasks/' model '/'];
    lambda    = task.lambda(1);
    mu        = task.lambda(2);
    mu_tsp    = task.lambda(3);
    suffix    = task.suffix;

    % check if selected
    if ~exist([modelpath 'slog' suffix '/lambda=' num2str(lambda) '_mu=' num2str(mu)])
        error(['there is no selected curves for lambda=', num2str(lambda)]);
    end

    disp(['start TSP: lambda=' num2str(lambda) ', mu=' num2str(mu)]);
    load([modelpath 'slog' suffix '/lambda=' num2str(lambda) '_mu=' num2str(mu) '/data.mat']);
    M = []; N = []; O= [];
    if task.type == 'm'
        cornerfile = ['../Models/' model '/corners.mat'];
        if exist(cornerfile)
            load(cornerfile);
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

    % get ready
    Amat(X==0, :) = [];
    Amat(:, X==0) = [];
    avail = reshape(Amat, [], 1);
    L2  = sum(X);
    nL2 = (L2 + 1) * (L2 + 1);
    k   = length(nums);

    %[f, A, b, Aeq, beq] = build_TSP(nums, E, D, lamda);   
    % [D E] = getcheat(D, E, 41, 1);

    if task.type == 's'
        [f, A, b, Aeq, beq, AvailIdx] = ...
            build_TSP2(k, E, D, mu_tsp, avail);  
        lb = [zeros(nL2, 1); ones(k, 1)];
        ub = [ones(nL2, 1);  ones(k, 1) * k];
        tic;
        x0 = gurobi_intlinprog(f, 1:size(f,1), A, b, Aeq, beq, lb, ub);  
        toc
    else
        %[f, A, b, Aeq, beq, quadcon, AvailIdx] = ...
        %    build_mTSP2(M, nums, E, D, mu_tsp, avail, task.mtsp_t);  
        maxCost = max(max(D + mu_tsp * E))
        
        for ik = 1: size(N, 2)
            E(N(1, ik), N(2, ik)) = 0;
            E(N(2, ik), N(1, ik)) = 0;
            D(N(1, ik), N(2, ik)) = 0;
            D(N(2, ik), N(1, ik)) = 0;
        end

        for ik = 1: size(M, 2)
            E(M(1, ik), M(2, ik)) = 100000;
            E(M(2, ik), M(1, ik)) = 100000;
            D(M(1, ik), M(2, ik)) = 100000;
            D(M(2, ik), M(1, ik)) = 100000;
        end

        if length(O)>0

            X(O) = 0;
        
            C = C(X==1);
            D = D(X==1, X==1);
            E = E(X==1, X==1);
            clusters = clusters(X==1);
            nums = [1];
            for it = 2: length(clusters)
                if clusters(it) == clusters(it-1)
                    nums(end) = nums(end) + 1;
                else
                    nums(end+1) = 1;
                end
            end
            k        = length(nums);
            L2       = sum(X);
            nL2      = (L2 + 1) * (L2 + 1);   
            avail    = ones(nL2, 1);
            Amat  = reshape(avail, L2+1, L2+1);
            F     = F(X==1); % new F
            X     = X(X==1); % new X
        end
        
        tic
        if given_m
            [f, A, b, Aeq, beq, AvailIdx] = ...
                build_mTSP(task.m, nums, E, D, mu_tsp, avail);  
            lb = [zeros(nL2, 1); ones(k, 1)];
            ub = [ones(nL2, 1);  ones(k, 1) * k];
        
            x0 = gurobi_intlinprog(f, 1:size(f,1), A, b, Aeq, beq, lb, ub);   
        else
            [f, A, b, Aeq, beq, quadcon, AvailIdx] = ...
                build_mTSP2(nums, E, D, mu_tsp, avail, 0.1 * maxCost);  
            lb = [zeros(nL2, 1); ones(k, 1); 1];
            ub = [ones(nL2, 1);  ones(k, 1) * k; 3];
            x0 = gurobi_qc(f, 1:size(f, 1), A, b, Aeq, beq, lb, ub, quadcon); 
        end
    end

    %% solve the integer linear programming using Gurobi interface.

    x=x0;
    %% prune some edges
    x=zeros(nL2+k,1);
    x(AvailIdx)=x0;


    % get orders and edges
    xij           = x(1: nL2);
    orders        = x(nL2+1: nL2+k);
    sxij          = find(xij == 1);
    edges         = [mod(sxij-1, k+1)+1, ceil(sxij / (k+1))];


    %% read-out selected points:
    source_path  = ['./tasks/' model '/new_off' suffix '/'];
    if ~exist([logdir '/x1'])
        mkdir([logdir '/x1']);
    end
    delete([logdir '/x1/*.off'])
    F1 = F(X==1);
    for it = 1: length(F1)   
          copyfile([source_path F1{it}], [logdir '/x1/' F1{it}]);
    end
    pointss = [];
    pointse = [];
    for it = 1: length(F1)
        points{it} = read_coff([logdir '/x1/' F1{it}]);     % read the points only once in the beginning.
        points{it} = points{it}(1:3, :);
        pointss = [pointss, points{it}(:, 1)];
        pointse = [pointse, points{it}(:, end)];
    end

    %% use curve distance
    D1 = pdist2(pointss', pointse');
    D2 = pdist2(pointse', pointss');
    D3 = pdist2(pointse', pointse');
    D4 = pdist2(pointss', pointss');

    % visualization here
    if verbose
        visualize(points,  1: length(F1), [], 'before TSP optimization', logdir);
        %visualize(points,  orders, [], 'after TSP optimization', logdir);
        visualize(points,  orders, edges, 'TSP edges', logdir);
    end

    disp('save TSP infomation')
    save([logdir '/tsp.mat'], 'edges', 'orders', 'points', 'D1', 'D2', 'D3', 'D4')


    %delete([logdir '/TSP_info.txt']);
    %fileID = fopen([logdir '/TSP_info.txt'],'w');
    %fprintf(fileID,'model\t %s\n',model);
    %fprintf(fileID,'lambda=%f, lamda=%f, tsp lambda=%f\n',lambda,lambda,mu_tsp);
    %fprintf(fileID,'TSP Time:\t %f sec\n',elapsedTime2);
    %fclose(fileID);

    alltime = toc
    disp('TSP done.')

    if task.cleantsp
        cleanTSP(task, false);
    end
end