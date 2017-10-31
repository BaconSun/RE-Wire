clc; clear;
warning('OFF');
addpath('./functions')

baseline_path = './line+tsp/'

% pika parameters
model    = 'fatcat'
mu       = 0.01;
lambda   = 1;
tsplam   = 0.015;
tol      = 5.e-7;
suffix   = '_allIter';


logdir = [baseline_path model]

% read points
[points lines] = read_obj([baseline_path model '/' model '-line3d-10.obj']);
%points([1, 2, 3], :) = points([3, 1, 2], :);


% curves
k = size(lines, 2);
nL2 = (k+1)^2;

for it = 1: k
    points1{it} = points(:, lines(:, it));
end


 % compute D & E
for x = 1: k
    for y = 1: k
        if x == y
            E(x, y) = 0; D(x, y) = 0; continue;
        end
        Dist    = pdist2(points1{x}', points1{y}');
        [r, c]  = find(Dist==min(min(Dist)));

        endDir1 = normalize(points1{x}(:, 3-r) - points1{x}(:, r));
        endDir2 = normalize(points1{y}(:, c) - points1{x}(:, 3-c));

        D(x, y) = Dist(r, c);
        E(x, y) = (1 - endDir1' * endDir2)/2;
        V(x, y, :) = [r, c];
    end
end

% build TSP
tic
[f, A, b, Aeq, beq, AvailIdx] = build_TSP2(k, E, D, tsplam, reshape(ones(k+1, k+1), [], 1)); 
lb = [zeros(nL2, 1); ones(k, 1)]; ub = [ones(nL2, 1);  ones(k, 1) * k];

%% prune some edges
f  = f(AvailIdx); A  = A(:,AvailIdx); Aeq= Aeq(:,AvailIdx); lb = lb(AvailIdx); ub = ub(AvailIdx);

% run TSP
x0 = gurobi_intlinprog(f, 1:size(f,1), A, b, Aeq, beq, lb, ub);  
elapsedTime2  = toc

x=zeros(nL2+k,1); x(AvailIdx)=x0;


% get orders and edges
xij           = x(1: nL2);
orders        = x(nL2+1: nL2+k);
sxij          = find(xij == 1);
edges         = [mod(sxij-1, k+1)+1, ceil(sxij / (k+1))];

[value, tour]= sort(orders);

% output starting/ending points (moving directions)
direction    = zeros(size(tour));
for it = 1: length(tour) - 1
    a = tour(it);
    b = tour(it + 1);
    v = V(a, b, :);
    if it == 1
        direction(it) = (v(1) == 1);
    end
    direction(it + 1) = (v(2) == 1);
end

% visualization here
visualize(points1,  1: k, [], 'before TSP optimization', logdir);
visualize(points1,  orders, [], 'after TSP optimization', logdir);
visualize(points1,  orders, edges, 'TSP edges', logdir);

points = points1;

disp('save TSP infomation')
save([logdir '/tsp.mat'], 'points', 'orders', 'direction', 'tour', 'logdir')

% fitting
% fit.
tic
[sp hp] = fitting(points, direction, tour, logdir, tol);
elapsedTime2 = toc

