function sd = cleaner(task, points1, pointss, pointse, sd, mark, use_floating, verbose)

    task.maxc = task.defloat(1:2);
    task.supporting = task.desupport(1:2);
    task.cutoff = task.defloat(3:4);
    task.d_th   = task.desupport(3:4);
    task.a_th   = task.desupport(5:6);
    task.s_th   = 0.01;
    task.rd_th  = 10;
    task.ra_th  = 0.75;
    task.end_th = 10;
    if strcmp(task.model, 'turtle_real')==1
        task.ra_th = 0.01;
        task.end_th = 40;
    end

    model   = task.model;
    markers = ['+', 'o', 'x',  '^', 'v', '>', '<', 's'];
    colours = ['r', 'g', 'b', 'm', 'c', 'y', 'k'];
    colors  = [255 0 0; 0 0 255; 0 255 255; 255 255 0; 0 255 0; 255 0 255]/255;
    if strcmp(mark, 'before')==1
        mm = 1;
    else
        mm = 2;
    end

    %% ---- check not supporting      
    D1 = pdist2(pointss', pointse') + eye(length(points1)) * 1000;
    D2 = pdist2(pointse', pointss') + eye(length(points1)) * 1000;
    D3 = pdist2(pointse', pointse') + eye(length(points1)) * 1000;
    D4 = pdist2(pointss', pointss') + eye(length(points1)) * 1000;

    for ik = 1: length(points1)
       curve = points1{ik};  
       [dd(ik, 1), dp(ik, 1)] = min(min(D1(ik, :), D4(ik, :)));
       [dd(ik, 2), dp(ik, 2)] = min(min(D2(ik, :), D3(ik, :)));
    end

    if verbose
        XF=figure; 
        set(XF, 'Position', [0, 0, 1600, 1200]);
    
        subplot(241)
        plot(abs(dd(:, 1)-dd(:, 2))); hold on;
        plot(ones(size(dd(:, 1))) * task.supporting(mm), 'k--')
        title([mark ':distance of supporting']);
    end

    for ik = 1: length(points1)
        supporting = min([norm(points1{ik}(1:3, 1) - points1{ik}(1:3, end)), task.supporting(mm)]);
        %ndd(ik) = (abs(dd(ik, 1)-dd(ik, 2))) > supporting || (dp(ik, 1) == dp(ik, 2));
        ndd(ik) = (abs(dd(ik, 1)-dd(ik, 2))) > supporting;
    end
    
    andd = cumsum(ndd) .* ndd;
    if verbose
        subplot(242);
        visualize2(points1, andd, [], ...
                  [mark ':check out supporting'], ... 
                  [], true, true, true);
    end

    %% --------------------------------
    camerapath = ['../Models/' model '/'];
    c1 = importdata([camerapath 'curvepoints1.txt']);
    c1 = reshape(c1(2: end), 5, []);
    c1 = c1(1:3, :);
    ce = [];
    le = max(c1(1, :));
    for e = 1: le
        c11 = find(c1(1, :)==e);
        V   = [c1(3, c11(1)); c1(2, c11(1))];
        ce  = [ce V];
        V   = [c1(3, c11(end)); c1(2, c11(end))];
        ce  = [ce V];
    end
    c_dist = squareform(pdist(ce'))+ eye(2 * le)*100000000;
    is_end_point = (min(c_dist) > task.end_th);

    camerapath = ['../Models/' model '/'];
    c2 = importdata([camerapath 'curvepoints2.txt']);
    c2 = reshape(c2(2: end), 5, []);
    c2 = c2(1:3, :);
    ce2 = [];
    le = max(c2(1, :));
    for e = 1: le
        c22  = find(c2(1, :)==e);
        V    = [c2(3, c22(1)); c2(2, c22(1))];
        ce2  = [ce2 V];
        V    = [c2(3, c22(end)); c2(2, c22(end))];
        ce2  = [ce2 V];
    end
    c_dist = squareform(pdist(ce2'))+ eye(2 * le)*100000000;
    is_end_point2 = (min(c_dist) > task.end_th);


    camerapath = ['../Models/' model '/'];
    c3= importdata([camerapath 'curvepoints3.txt']);
    c3 = reshape(c3(2: end), 5, []);
    c3 = c3(1:3, :);
    ce3 = [];
    le = max(c3(1, :));
    for e = 1: le
        c33  = find(c3(1, :)==e);
        V    = [c3(3, c33(1)); c3(2, c33(1))];
        ce3  = [ce3 V];
        V    = [c3(3, c33(end)); c3(2, c33(end))];
        ce3  = [ce3 V];
    end
    c_dist = squareform(pdist(ce3'))+ eye(2 * le)*100000000;
    is_end_point3 = (min(c_dist) > task.end_th);


    T1 = read_curvepoints([camerapath 'curvepoints1.txt']);
    try
        M1 = read_camera([camerapath 'camera1.txt']);   
    catch
        M1 = read_camera2([camerapath 'camera1.txt']);
    end
    T2 = read_curvepoints([camerapath 'curvepoints2.txt']);
    try
        M2 = read_camera([camerapath 'camera2.txt']);   
    catch
        M2 = read_camera2([camerapath 'camera2.txt']);
    end
    T3 = read_curvepoints([camerapath 'curvepoints3.txt']);
    try
        M3 = read_camera([camerapath 'camera3.txt']);   
    catch
        M3 = read_camera2([camerapath 'camera3.txt']);
    end


    if verbose
        subplot(243);        
        
        %% --distanec
        x = -1:0.01:1;
        y = task.d_th(mm) * ones(size(x));
        plot(y, x, 'k', 'linewidth', 2); hold on;
        
        x = 0:0.01:1;
        y = task.a_th(mm) * ones(size(x));
        plot(x, y, 'k', 'linewidth', 2); hold on;
    end

    num = 0;
    for jk = 1: length(points1)
        if ndd(jk) == 0
            u(jk) = 1;
            continue;
        end
        num = num+1;
    
        %% -- angle
        L  = min(size(points1{jk}, 2), 5);
        
        if dd(jk, 1) > dd(jk, 2)  % start not is supporting
            Vj = points1{jk}(1:3, 1)   - points1{jk}(1:3, L);
        else                      % end is not supporting
            Vj = points1{jk}(1:3, end) - points1{jk}(1:3, end-L+1);
        end
        
        for ik = 1: length(points1)
            if dd(jk, 1) > dd(jk, 2)  % start is not supporting
                Uj(1, :) = points1{ik}(1:3, 1)   - points1{jk}(1:3, 1);
                Uj(2, :) = points1{ik}(1:3, end) - points1{jk}(1:3, 1);
                Dj(1, :) = D4(jk,:); % s-e
                Dj(2, :) = D1(jk,:); % s-s
            else
                Uj(1, :) = points1{ik}(1:3, 1)   - points1{jk}(1:3, end);
                Uj(2, :) = points1{ik}(1:3, end) - points1{jk}(1:3, end);
                Dj(1, :) = D2(jk,:); % e-s
                Dj(2, :) = D3(jk,:); % e-e
            end
        
            Aj(1, ik) = dot(Vj, Uj(1, :))/(norm(Vj, 2) * norm(Uj(1, :), 2));
            Aj(2, ik) = dot(Vj, Uj(2, :))/(norm(Vj, 2) * norm(Uj(2, :), 2));
        end
        
        support = sum((Dj < task.d_th(mm)) & (Aj > task.a_th(mm)), 2)>0;
        indexs  = find(Dj < 1000);

        
        if verbose
            for sk=1:length(indexs)
                %if (Dj(indexs(sk)) < task.d_th(mm) + 0.2) & (Aj(indexs(sk)) > task.a_th(mm)-0.2)
                %if (num==3)   
                if strcmp(mark, 'before')~=1
                    text(Dj(indexs(sk)),Aj(indexs(sk)), num2str(num));hold on;
                end
                %end
            end
            xlim([0, 0.2]); ylim([0, 1])
            xlabel('distance');ylabel('cos(theta)');
        end

        % if (num == 3)
        %     disp('--------')
        %     jk
        %     support
        % end

        % ---- compute the angle -----
        if (support(1)) || (support(2))
            u(jk)=1;
        else
            u(jk)=0;


            if (num == 3)
                u(jk)
            end

            % project to the image and look for the end-point      
            if dd(jk, 1) > dd(jk, 2) % start point farther (not supporting)
                proj = points1{jk}(1:3, 1);
            else                     % end point farther
                proj = points1{jk}(1:3, end);
            end
            U = project(M1, proj);
            dis_to_end = sqrt(sum((ce - repmat(U', [1, size(ce, 2)])).^2, 1));
            [min_dis, pos] = min(dis_to_end);
            if is_end_point(pos) == 1
                if min_dis < 30
                    u(jk)=1;
                end
            end

            % if (num == 3)
            %     dis_to_end
            %     u(jk)
            % end

            U = project(M2, proj);
            dis_to_end = sqrt(sum((ce2 - repmat(U', [1, size(ce2, 2)])).^2, 1));
            [min_dis, pos] = min(dis_to_end);
            if is_end_point2(pos) == 1
                if min_dis < 30
                    u(jk)=1;
                end
            end

            % if (num == 3)
            %     dis_to_end
            %     u(jk)
            % end

            U = project(M3, proj);
            dis_to_end = sqrt(sum((ce3 - repmat(U', [1, size(ce3, 2)])).^2, 1));
            [min_dis, pos] = min(dis_to_end);
            if is_end_point3(pos) == 1
                if min_dis < 30
                    u(jk)=1;
                end
            end

            % if (num == 3)
            %     dis_to_end
            %     u(jk)
            % end

        end
    end  

    if verbose
        subplot(244);
        visualize2(points1, u, [], ...
                  [mark ':delete supporting'], ... 
                  [], false, colors, true);
    end

    selected_edges = u';

    %% -- recompute data (copy from preprocess)     
    X1_ = find(sd.X == 1);
    X1_ = X1_(selected_edges==1);
    sd.X   = sd.X * 0;
    sd.X(X1_) = 1;   %% new X
    
    a = 0;
    points_ = {}; pointss1 = []; pointse1 = [];
    for ix = 1: length(sd.X)
        if sd.X(ix) == 1
            a = a + 1;
            points_{a} = points1{ix};
            pointss1 = [pointss1, points_{a}(1:3, 1)];
            pointse1 = [pointse1, points_{a}(1:3, end)];
        end
    end

    points1 = points_;
    sd.points  = points1;
    sd.pointss = pointss1;
    sd.pointse = pointse1;
    
    sd.C = sd.C(sd.X==1);
    sd.D = sd.D(selected_edges==1, selected_edges==1);
    sd.E = sd.E(selected_edges==1, selected_edges==1);
    sd.clusters = sd.clusters(sd.X==1);
    nums = [1];
    for it = 2: length(sd.clusters)
        if sd.clusters(it) == sd.clusters(it-1)
            nums(end) = nums(end) + 1;
        else
            nums(end+1) = 1;
        end
    end

    k        = length(nums);
    L        = sum(nums);
    nL       = (L + 1) * (L + 1);   
    avail    = ones(nL, 1);
    sd.Amat  = reshape(avail, L+1, L+1);
    sd.F     = sd.F(sd.X==1); % new F
    sd.X     = sd.X(sd.X==1); % new X
    sd.nums  = nums;


    if use_floating
        %% -- delete floating curves --
        G = sd.D;
        Z = linkage(G, 'single');
        % dendrogram(Z)
        T = cluster(Z, 'maxclust', task.maxc(mm));

        numc = max(T);
        Gc = zeros(numc, numc);
        G_row = [];
        G_col = [];
        for t1 = 1: numc
            for t2 = 1: numc
                if t1 == t2
                    continue
                end
                Gc(t1, t2) = min(min(G(T==t1, T==t2)));
            end
        end
        
        for k = 1: numc
            J(k) = min(Gc(k, (1:numc) ~= k));
        end

        if verbose
            subplot(245)
            plot(J); hold on;
            plot(ones(size(J)) * task.cutoff(mm), 'k--'); hold on;
            title([mark ':distance of floating clusters']);

            subplot(246);
            visualize2(points1, T, [], ...
                      [mark ':clustering results'], ... 
                      [], false, true, true);
        end
        selected_clusters = find(J < task.cutoff(mm));
        selected_edges = ismember(T, selected_clusters);
        
        if verbose    
            subplot(247);
            visualize2(points1, selected_edges, [], ...
                      [mark ':select floating curves'], ... 
                      [], false, colors, true);

            %% -- merge
            %% selected_edges = selected_edges .* u';   
            subplot(248);
            visualize2(points1, selected_edges, [], ...
                      [mark ':all removed'], ... 
                      [], false, colors, true);
        end
        
        %% -- recompute data (copy from preprocess)     
        X1_ = find(sd.X == 1);
        X1_ = X1_(selected_edges==1);
        sd.X   = sd.X * 0;
        sd.X(X1_) = 1;   %% new X

        sd.C = sd.C(sd.X==1);
        sd.D = sd.D(selected_edges==1, selected_edges==1);
        sd.E = sd.E(selected_edges==1, selected_edges==1);
        sd.clusters = sd.clusters(sd.X==1);
        nums = [1];
        for it = 2: length(sd.clusters)
            if sd.clusters(it) == sd.clusters(it-1)
                nums(end) = nums(end) + 1;
            else
                nums(end+1) = 1;
            end
        end

        a = 0;
        points_ = {}; pointss1 = []; pointse1 = [];
        for ix = 1: length(sd.X)
            if sd.X(ix) == 1
                a = a + 1;
                points_{a} = points1{ix};
                pointss1 = [pointss1, points_{a}(1:3, 1)];
                pointse1 = [pointse1, points_{a}(1:3, end)];
            end
        end
        points1    = points_;
        sd.points  = points1;
        sd.pointss = pointss1;
        sd.pointse = pointse1;

        k        = length(nums);
        L        = sum(nums);
        nL       = (L + 1) * (L + 1);   
        avail    = ones(nL, 1);
        sd.Amat  = reshape(avail, L+1, L+1);
        sd.F     = sd.F(sd.X==1); % new F
        sd.X     = sd.X(sd.X==1); % new X
        sd.nums  = nums;
    end
end