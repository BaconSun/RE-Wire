function finalfit(task, verbose)

    baseline_path = './line+tsp/';
    model  = task.model;
    modelname = model;
    suffix = task.suffix;
    tol    = task.tol(1);
    ctol   = task.tol(2);
    lambda = task.lambda(1);
    mu     = task.lambda(2);

    if strcmp(task.model,'turtle')==1   % special setting for turtle
        tolx  = [3.e+6, 1.5 *1.e+5, 3.e+6];
    else
        tolx  = [ctol ctol ctol];
    end

    modelpath = ['./tasks/' model '/'];
    logdir    = [modelpath 'log' ... 
                 '/lambda=' num2str(lambda) '_mu=' num2str(mu)];
    showIter  = true;
    iters     = 15;

    load([logdir '/tsp.mat'])
    %load([baseline_path 'fatcat/tsp.mat'])

    % points check
    for it = 1: length(points)
        points{it} = points{it}(1:3, :);
    end
    
    finalFig = figure;

    % get tours from edges
    if verbose    
        set(finalFig, 'Position', [100, 100, 800, 600]);
    end

    O      = max(max(edges));
    starts = find(edges(:, 1) == O);
    fp_all = [];
    alphas = 0: 180; alphas = alphas/180;
    tic
    for s = 1:  length(starts)
        tour = [];
        s1   = starts(s);
        while true
            v1   = edges(s1, 2);
            tour = [tour v1];
            s1 = find(edges(:, 1) == v1);
            if edges(s1, 2) == O
                break
            end
        end
        
        % output starting/ending points (moving directions)
        direction    = zeros(size(tour));
        for it = 1: length(tour) - 1
            a = tour(it);
            b = tour(it + 1);
            d1 = D1(a, b); d2 = D2(a, b); d3 = D3(a, b); d4 = D4(a, b);
            if it == 1
                direction(it) = (min(d4, d1) > min(d2, d3));
            end
            direction(it + 1) = (min(d1, d3) > min(d2, d4));
        end


        

        for k=[1 length(tour)]
            ts = tour(k);
            if direction(k) == 1 % forward direction
                ds = points{ts}(:, 1);
            else
                ds = points{ts}(:, end);
            end
            
            if exist('T')
                % connect lines.
                st = find(T(:, 1) == ts)
                if st > 0
                    fs = T(st, 2);
                    d1 = norm(ds - points{fs}(:, 1));
                    d2 = norm(ds - points{fs}(:, end));

                    if d1 < d2
                        ps = points{fs}(:, 1);
                    else
                        ps = points{fs}(:, end);
                    end
                    if ~xor(direction(k) == 1, k==1) % forward
                        p  = ds * (alphas - 0.17) + ps * (1 - alphas + 0.17);
                        points{ts} = [p points{ts}];
                    else
                        p  = ds * (alphas + 0.17) + ps * (1 - alphas - 0.17);
                        points{ts} = [points{ts} fliplr(p)];
                    end
                end
            end
        end
        
        % fitting
        [sp{s} fp{s}]= fitting(finalFig, points, direction, tour, logdir, tol, verbose);
        fp_all = [fp_all fp{s}];
    end


    % ------------------------------------------------------------
    % load([baseline_path 'fatcat/tsp.mat'])
    % FF = figure;
    % [sp hp] = fitting(points, direction, tour, logdir, tol, FF);
    % return
    % final fitting with images
    % read image points!
    % ------------------------------------------------------------

    final_points = [];
    group_points = {};
    colors = ['r', 'b', 'g', 'y', 'm'];
    camerapath = ['../Models/' modelname '/'];
    for s = 1:  length(starts)
        pp   = size(fp{s}, 2) / size(fp_all, 2);
        pnum = 400000 * pp;
        [result{s} dd{s}] = togetherfit(camerapath, sp{s}, tol, pp * tolx(s), pnum, 10, verbose);
        group_points{s} = result{s};
        final_points = [final_points result{s}];
    end
    alltime = toc

    figure;
    for s = 1: length(starts)
        plot3(group_points{s}(1,:), group_points{s}(2,:), group_points{s}(3,:),...
              ['.' colors(s)]); hold on;
    end
    axis equal;
    hold off;
    
    outdir = ['results/' model];
    write_off([outdir '/final_' model(1:end-1) '.off'], final_points, int16.empty(0,3));
    for s = 1: length(starts)
        write_off([outdir '/final_' model(1:end-1) '_curve_' num2str(s) '.off'], result{s}, int16.empty(0,3));
    end
    


    % save the final curve
    %delete([logdir '/Fitting_info.txt']);
    %fileID = fopen([logdir '/Fitting_info.txt'],'w');
    %fprintf(fileID,'model\t %s\n',model);
    %fprintf(fileID,'tol=1/%f, ctol=%f\n',1/tol,ctol);
    %fprintf(fileID,'Fitting 3D PointsTime:\t %f sec\n' ,elapsedTime2);
    %fprintf(fileID,'Fitting Views Time:\t %f sec\n',elapsedTime3);
    %fclose(fileID);
    %disp('Save Fitting information');
end


