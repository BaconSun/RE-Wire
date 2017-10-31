function preprocess(task, verbose)
    recompute_DE = false;
    delete_node  = true;

    model  = task.model;
    suffix = task.suffix;
    modelpath = ['./tasks/' model '/'];

    % before doing everything:
    source_path  = ['../Models/' model '/notiso_color'  suffix '/'];
    F = importdata([modelpath 'CurveFileName' suffix '.txt']);
    W = importdata([modelpath 'outCEachView.txt']);
    for it = 1: length(F)   
          points0{it} = read_coff([source_path F{it}]);
          curveNums0(it) = size(points0{it}, 2);
          points0{it}(3, :) = -points0{it}(3, :);
          %points0{it}(2, :) = -points0{it}(2, :);
    end

    % get ready
    nums = importdata([modelpath 'outPartNum' suffix '.txt']);   % input the number of each parts
    nums = nums';
    nums = nums(nums>0);
    k        = length(nums);
    L        = sum(nums);
    nL       = (L + 1) * (L + 1);
    cc       = [0 cumsum(nums)];                                         % cluster boundaries
    clusters = [];
    for n = 1: k
        clusters = [clusters; ones(nums(n), 1) * n];                     % node-id  -> cluster-id
    end
    avail = ones(nL, 1);
    Amat = reshape(avail, L+1, L+1);

    % load C/D/E
    C    = importdata([modelpath 'outC' suffix '.txt']);
    D    = importdata([modelpath 'outD' suffix '.txt']);
    E    = importdata([modelpath 'outE' suffix '.txt']);
    if recompute_DE
        [E, D] = recompute(source_path, F);
    end

    % visualization
    if verbose
        [H, colors, fH, T]   = visualize3(points0, clusters, C, F, [], 'before selection', ... 
                                   modelpath, false, true); % all the curves
        axis equal; view(-38, -20);
        %set(H, 'ButtonDownFcn', {@selected2, H, clusters, T});
        %waitfor(fH);
    end
    X    = ones(length(F), 1);

    if delete_node
        oldIDs     = find(X==1);
       
        if ~exist([modelpath 'new_off' suffix])
            mkdir([modelpath 'new_off' suffix]);
        end

        F = F(X==1);
        for it = 1: length(F)   
            copyfile([source_path F{it}], [modelpath 'new_off' suffix '/' F{it}]);
            points{it} = points0{oldIDs(it)};
        end
        source_path  = [modelpath 'new_off' suffix '/'];


        % recompute new data
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
        L        = sum(nums);
        nL       = (L + 1) * (L + 1);   
        avail    = ones(nL, 1);
        Amat     = reshape(avail, L+1, L+1);
    end

    save([modelpath 'pre_data' suffix '_.mat'], 'C', 'D', 'E', 'F', ... 
                    'Amat', 'nums', 'clusters', 'source_path', 'suffix');
    %save(['./cheat/' task.cheatfile], 'deletedIDs');
    disp('done.')

end