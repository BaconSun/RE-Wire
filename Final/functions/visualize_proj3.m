function [TS0, TS, TT, VS] = visualize_proj3(t1, M1, points, points2, orders, features, edges, ...
                                             name, logdir, showtext, showColor, model)
    rand_colors = rand(500, 3);
    % if not given orders
    if length(orders) == 0
        orders = 1: size(points, 2);
    end

    FigHandle = figure;
    set(FigHandle, 'Position', [0, 0, 1000, 2000]);
    pointsm = [];
 
    
    % subplot(211); 
    TT = [];
    TD = [];
    TS = [];    
    for pi =1: max(t1(1, :))
        T1{pi} = [];
    end
    for i = 1: size(t1, 2)
        T1{t1(1, i)} = [T1{t1(1, i)} t1(3:-1:2, i)];
        TT = [TT t1(3:-1:2, i)];
        TD = [TD t1(1, i)];
        TS = [TS 0];
    end
    KNNTT = createns(TT', 'nsmethod', 'kdtree');

    for it = 1: length(points)   
        V1 = project(M1, points{it});
        VD{it} = [];
        
        % find the neraest curve point
        [idx dist] = knnsearch(KNNTT, V1);
        curve_ids  = TD(idx);
        unique_ids = unique(curve_ids);
        VS{it} = [];
        for jt = 1: length(unique_ids)
            uid = unique_ids(jt);
            rid = idx(curve_ids == uid);
            min_rid = min(rid);
            max_rid = max(rid);
            TS(min_rid: max_rid) = TS(min_rid: max_rid) + 1;
            % plot(TT(1, min_rid: max_rid), TT(2, min_rid: max_rid)); hold on;
            VS{it} = [VS{it}; [uid, min_rid, max_rid]];
        end   
    end
    
    % task 1: show all double cover
    covers = [];
    Vdt = {};
    for ik = 1: length(TD)
        cover = [];
        for iu = 1: length(VS)
            for iq = 1: size(VS{iu}, 1)
                if ((ik > (VS{iu}(iq, 2)-1)) && (ik < (VS{iu}(iq, 3)+1)))
                    cover = [cover; iu];
                end
            end
        end
        for iv = 1: length(cover)
            VD{cover(iv)} = [VD{cover(iv)}; cover];
        end
        Vdt{ik} = cover;
        
        if length(cover) > 1    
            covers = [covers; cover];
        end
    end
    
    for id = 1: length(VD)
        VD{id} = unique(VD{id});
    end
    plot(TT(1, :), TT(2, :), '.r'); hold on;
    %plot(TT(1, :), TT(2, :), '.r'); hold on;
    plot(TT(1, TS>0), TT(2, TS>0), '.b'); hold on;
    axis equal;
    
    TS0 = TS;
   
    TT = [];
    TD = [];
    TS = [];    
    for pi =1: max(t1(1, :))
        T1{pi} = [];
    end
    for i = 1: size(t1, 2)
        T1{t1(1, i)} = [T1{t1(1, i)} t1(3:-1:2, i)];
        TT = [TT t1(3:-1:2, i)];
        TD = [TD t1(1, i)];
        TS = [TS 0];
    end
    KNNTT = createns(TT', 'nsmethod', 'kdtree');

    VS = {};
    VD = {};
    
    for it = 1: length(points2)   
        V1 = project(M1, points2{it});
        VD{it} = [];
        
        % find the neraest curve point
        [idx dist] = knnsearch(KNNTT, V1);
        curve_ids  = TD(idx);
        unique_ids = unique(curve_ids);
        VS{it} = [];
        for jt = 1: length(unique_ids)
            uid = unique_ids(jt);
            rid = idx(curve_ids == uid);
            min_rid = min(rid);
            max_rid = max(rid);
            TS(min_rid: max_rid) = TS(min_rid: max_rid) + 1;
            % plot(TT(1, min_rid: max_rid), TT(2, min_rid: max_rid)); hold on;
            VS{it} = [VS{it}; [uid, min_rid, max_rid]];
        end   
    end
    
    % task 1: show all double cover
    covers = [];
    Vdt = {};
    for ik = 1: length(TD)
        cover = [];
        for iu = 1: length(VS)
            for iq = 1: size(VS{iu}, 1)
                if ((ik > (VS{iu}(iq, 2)-1)) && (ik < (VS{iu}(iq, 3)+1)))
                    cover = [cover; iu];
                end
            end
        end
        for iv = 1: length(cover)
            VD{cover(iv)} = [VD{cover(iv)}; cover];
        end
        Vdt{ik} = cover;
        
        if length(cover) > 1    
            covers = [covers; cover];
        end
    end
    
    for id = 1: length(VD)
        VD{id} = unique(VD{id});
    end
    
    %view(-6, -34)
    %savefig([logdir '/' name '.fig']);
end

