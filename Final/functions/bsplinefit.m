modelname='turtle_real'

addpath('functions/');
addpath('functions/bspline/');
baseline_path = './line+tsp/';

task   = setup(modelname);
%task   = setup('cart');
model  = task.model;
suffix = task.suffix;
tol    = task.tol;
ctol   = task.ctol;
mu     = task.select(1);
l      = task.select(2);

modelpath = ['./tasks/' model '/'];
logdir    = [modelpath 'log' ... 
             '/mu=' num2str(mu) '_l=' num2str(l)];
showIter  = true;
iters     = 15;

load([logdir '/tsp.mat'])
%load([baseline_path 'fatcat/tsp.mat'])

% points check
for it = 1: length(points)
    points{it} = points{it}(1:3, :);
end


% get tours from edges
finalFig = figure;
set(finalFig, 'Position', [100, 100, 800, 600]);

O      = max(max(edges));
starts = find(edges(:, 1) == O);
fp_all = [];
savepoints = [];

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
    [sp{s} fp{s} pp{s}]= fitting(finalFig, points, direction, tour, logdir, tol);
    savepoints = [savepoints; pp{s}];
    fp_all = [fp_all fp{s}];
end

write_off([logdir '/before_image_fit.off'], savepoints, int16.empty(0,3));
