function filter_curves(task, options, verbose)
tic;
ifcheck     = false;
beforeQO    = true;
afterQO     = true;
checkcover  = false;
min_overlap = 5;

lambda     = task.lambda(1);
mu         = task.lambda(2);
suffix     = task.suffix;
model      = task.model;
modelpath = ['./tasks/' model '/'];

% load merged data
load([modelpath 'pre_data' suffix '_.mat']);
suffix = task.suffix;
if ~exist([modelpath 'slog' suffix])
    mkdir([modelpath 'slog' suffix]);
end
camerapath = ['../Models/' model '/'];
source_path  = ['../Models/' model '/notiso_color' suffix '/']
for it = 1: length(F)   
      points0{it} = read_coff([source_path F{it}]);
end

length(points0);

if ifcheck
    set(H, 'ButtonDownFcn', {@check_cluster, H});
    waitfor(fH);
end


% log-files
logdir=[modelpath 'slog' suffix '/lambda=' num2str(lambda) '_mu=' num2str(mu)];
if ~exist(logdir)
    mkdir(logdir);
end
X  = ones(sum(nums), 1);

sd = {};
sd.X = X; sd.F=F; sd.C=C;
sd.D = D; sd.E=E; sd.Amat = Amat;
sd.clusters = clusters;
pointss = []; pointse = [];
for it = 1: length(F)   
      points0{it} = read_coff([source_path F{it}]);
      % points0{it}(3, :) = -points0{it}(3, :); 
      pointss = [pointss, points0{it}(1:3, 1)];
      pointse = [pointse, points0{it}(1:3, end)];
end


if beforeQO
    sd2 = sd;
    remained = length(points0)
    while true
        sd2 = cleaner(task, points0, pointss, pointse, sd2, 'before', false, verbose);
        X = sd2.X; F=sd2.F; C=sd2.C; D = sd2.D; E=sd2.E; 
        Amat = sd2.Amat; clusters = sd2.clusters; nums = sd2.nums;
        points0 = sd2.points; pointss = sd2.pointss; pointse = sd2.pointse;
        if length(points0) >= remained
            break;
        end
        remained = length(points0)
    end
    
    % with floating
    sd2 = cleaner(task, points0, pointss, pointse, sd2, 'before', true, verbose);
    X = sd2.X; F=sd2.F; C=sd2.C; D = sd2.D; E=sd2.E; 
    Amat = sd2.Amat; clusters = sd2.clusters; nums = sd2.nums;

end

% Round-up
[X, E, D, O]=Roundup(C, E, D, nums, lambda, mu);

%% read-out selected points:
if ~exist([logdir '/x1' suffix])
    mkdir([logdir '/x1' suffix]);
end
delete([logdir '/x1' suffix '/*.off'])
F = F(X==1);
clusters = clusters(X==1);
points  = [];
pointss = [];
pointse = [];
for it = 1: length(F)   
      copyfile([source_path F{it}], [logdir '/x1' suffix '/' F{it}]);
      points1{it} = read_coff([logdir '/x1' suffix '/' F{it}]);
      % points1{it}(3, :) = -points1{it}(3, :);
      points = [points points1{it}]; 
      pointss = [pointss, points1{it}(1:3, 1)];
      pointse = [pointse, points1{it}(1:3, end)];
end
length(points1)
if verbose
    disp(['save lambda = ' num2str(lambda) '_mu=' num2str(mu)]);
    visualize2(points1, 1:length(X==1), [], ...
              ['After QO' '-lambda=' num2str(lambda) '-mu=' num2str(mu)], ... 
              logdir, false, true);
    view(-33, -2); axis equal;
    view(-61, 40);
end

sd = {};
sd.X = X(X==1); sd.F=F; sd.C=C;
sd.D = D; sd.E=E; sd.Amat = Amat;
sd.clusters = clusters;
%return

if afterQO

    sd2 = sd;
    remained = length(points1)
    while true
        sd2 = cleaner(task, points1, pointss, pointse, sd2, 'after', false, verbose);
        X = sd2.X; F=sd2.F; C=sd2.C; D = sd2.D; E=sd2.E; 
        Amat = sd2.Amat; clusters = sd2.clusters; nums = sd2.nums;
        points1 = sd2.points; pointss = sd2.pointss; pointse = sd2.pointse;
        if length(points1) >= remained
            break;
        end
        remained = length(points1)
    end
    
    % with floating
    sd2 = cleaner(task, points1, pointss, pointse, sd2, 'after', true, verbose);
    X = sd2.X; F=sd2.F; C=sd2.C; D = sd2.D; E=sd2.E; 
    Amat = sd2.Amat; clusters = sd2.clusters; nums = sd2.nums;
end       

%% read-out selected points (again):
if ~exist([logdir '/x1' suffix])
    mkdir([logdir '/x1' suffix]);
end
delete([logdir '/x1' suffix '/*.off']);

points  = [];
points2 = {};
for it = 1: length(F)   
      copyfile([source_path F{it}], [logdir '/x1' suffix '/' F{it}]);
      points2{it} = read_coff([logdir '/x1' suffix '/' F{it}]);
      points2{it}(3, :) = -points2{it}(3, :);
      points = [points points2{it}]; 
end

if verbose
    disp(['save lambda = ' num2str(lambda) '_mu=' num2str(mu)]);
    visualize2(points2, 1:length(X==1), [], ...
              ['After Deleting All Curves' '-lambda=' num2str(lambda) '-mu=' num2str(mu)], ... 
              logdir, false, true);
    view(-33, -2); axis equal;
    view(-61, 40);
end

% return
% save a merged file, easy for read
write_coff([logdir '/selected_lambda=' num2str(lambda) '_mu=' num2str(mu) '.off'], points, int16.empty(0,3));
save([logdir '/data.mat'], 'X', 'C', 'D', 'E', 'F', 'Amat', ...
                           'nums', 'clusters', 'source_path');


alltime = toc
end