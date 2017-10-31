clc; clear; tic;
warning('OFF');
addpath('./functions')

ifcheck = false;
%task = setup('human_real');
%task = setup('horse_real5');
%task = setup('flower_real');
%task   = setup('bird_real');
%task   = setup('elephant_real');
%task   = setup('bike_real7');
%task   = setup('fatcat_real');
task   = setup('turtle_real');
mus    = task.mu;
lams   = task.lambda;
suffix = task.suffix;
model  = task.model;
modelpath = ['./tasks/' model '/'];

% load merged data
load([modelpath 'pre_data' suffix '_.mat']);
C0 = C; D0 = D; E0 =E; 
suffix = task.suffix;
if ~exist([modelpath 'slog' suffix])
    mkdir([modelpath 'slog' suffix]);
end
camerapath = ['../Models/' model '/'];
source_path  = ['../Models/' model '/notiso_color' suffix '/']
for it = 1: length(F)   
      points0{it} = read_coff([source_path F{it}]);
end

length(points0)

if ifcheck
    set(H, 'ButtonDownFcn', {@check_cluster, H});
    waitfor(fH);
end

for mu = mus(1): mus(2): mus(3)
    for l = lams(1): lams(2): lams(3)
        % log-files
        logdir=[modelpath 'slog' suffix '/mu=' num2str(mu) '_l=' num2str(l)];
        if ~exist(logdir)
            mkdir(logdir);
        end

        % gurobi QP
        % tic; [X, E, D]=gurobi_intqp(C0 , E0, D0, nums, mu, l);

        % Round-up
        [X, E, D, O]=Roundup(C0, E0, D0, nums, mu, l);
        
        %% read-out selected points:
        if ~exist([logdir '/x1' suffix])
            mkdir([logdir '/x1' suffix]);
        end
        delete([logdir '/x1' suffix '/*.off'])
        
        F1 = F(X==1);
        points  = [];
        for it = 1: length(F1)   
              copyfile([source_path F1{it}], [logdir '/x1' suffix '/' F1{it}]);
              points1{it} = read_coff([logdir '/x1' suffix '/' F1{it}]);
              points = [points points1{it}]; 
        end
        length(points1);
        disp(['save mu = ' num2str(mu) '_l=' num2str(l)]);
        [H, rndc, fH, I] = visualize_proj(points1, 1:length(X==1), [], [], ...
                  ['After Round-up, Nothing deleted' '-mu=' num2str(mu) '-l=' num2str(l)], ... 
                  logdir, false, true, camerapath);
        set(H, 'ButtonDownFcn', {@selected3, H, I});      
        waitfor(fH);
        return
        
        %% Just for visualization?
        
        
        %% -- delete curves --
        Z = linkage(D, 'single');
        dendrogram(Z)
        T = cluster(Z, 'maxclust', 8)
        visualize2(points1, T, [], ...
                  ['After Clustering ' '-mu=' num2str(mu) '-l=' num2str(l)], ... 
                  logdir, false, true);
        numc = max(T);
        Dc = zeros(numc, numc);
        for t1 = 1: numc
            for t2 = 1: numc
                Dc(t1, t2) = min(min(D(T==t1, T==t2)));
            end
        end
        for k = 1: numc
            J(k) = min(Dc(k, (1:numc) ~= k));
        end
        figure; plot(J)
        selected_clusters = find(J < 0.01);
        selected_edges = ismember(T, selected_clusters);
        visualize2(points1, selected_edges, [], ...
                  ['After Round-up ' '-mu=' num2str(mu) '-l=' num2str(l)], ... 
                  logdir, false, true);
        
        %% -- recompute data (copy from preprocess)     
        X1_ = find(X == 1);
        X1_ = X1_(selected_edges==1);
        X      = X * 0;
        X(X1_) = 1;   %% new X
        
        C = C(X==1);
        D = D(selected_edges==1, selected_edges==1);
        E = E(selected_edges==1, selected_edges==1);
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
        F        = F(X==1); % new F
        X        = X(X==1); % new X
        
        %% read-out selected points (again):
        if ~exist([logdir '/x1' suffix])
            mkdir([logdir '/x1' suffix]);
        end
        delete([logdir '/x1' suffix '/*.off'])
        
        F1 = F(X==1);
        points  = [];
        points2 = {};
        for it = 1: length(F1)   
              copyfile([source_path F1{it}], [logdir '/x1' suffix '/' F1{it}]);
              points2{it} = read_coff([logdir '/x1' suffix '/' F1{it}]);
              points = [points points2{it}]; 
        end
        length(points2);
        figure;
        disp(['save mu = ' num2str(mu) '_l=' num2str(l)]);
        visualize2(points2, 1:length(X==1), [], ...
                  ['After Deleting Curves' '-mu=' num2str(mu) '-l=' num2str(l)], ... 
                  logdir, false, true);
        view(-33, -2); axis equal;
        view(-61, 40);

        % save a merged file, easy for read
        write_coff([logdir '/selected_mu=' num2str(mu) '_l=' num2str(l) '.off'], points, int16.empty(0,3));
        save([logdir '/data.mat'], 'X', 'C', 'D', 'E', 'F', 'Amat', ...
                                   'nums', 'clusters', 'source_path');
    end
end

alltime = toc