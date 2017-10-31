function pick_corner(task)

    mu = task.lambda(1);
    l  = task.lambda(2);
    model = task.model;
    modelpath = ['./tasks/' model '/'];

    M = []; N = []; O = [];
    logdir=[modelpath 'slog/mu=' num2str(mu) '_l=' num2str(l)];
    load([logdir '/data.mat']);

    if exist([logdir '/corners.mat'])
        load([logdir '/corners.mat']); 
        M 
        N
        O
    end

    F1 = F(X==1);
    for it = 1: length(F1)   
          points1{it} = read_coff([logdir '/x1/' F1{it}]);
    end
    [H c fH] = visualize2(points1, 1: length(points1), [], ...
                          'Round-up', logdir, true, true);

    S = {}; % picked corners
    set(H, 'ButtonDownFcn', {@selected0, H});
    waitfor(fH);

    %M = [];
    for it = 1: length(S)
        if length(S{it}) > 0
            for k = 1: length(S{it})
                M = [M [it; S{it}(k)]];
            end
        end
    end
    %save([logdir '/corners.mat'], 'M')



    
    [H c fH] = visualize2(points1, 1: length(points1), [], ...
                          'Round-up', logdir, true, true);

    S = {}; % picked corners
    set(H, 'ButtonDownFcn', {@selected0, H});
    waitfor(fH);

    %M = [];
    for it = 1: length(S)
        if length(S{it}) > 0
            for k = 1: length(S{it})
                N = [N [it; S{it}(k)]];
            end
        end
    end

    [H c fH] = visualize2(points1, 1: length(points1), [], ...
                          'Round-up', logdir, true, true);

    X = []; % picked corners
    set(H, 'ButtonDownFcn', {@selected, H});
    waitfor(fH);

    %M = [];
    O = find(X==0);

    save([logdir '/corners.mat'], 'M', 'N', 'O')
end



