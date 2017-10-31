clc; clear;
warning('OFF');
addpath('./functions/');

model  = 'cart/'
mu     = 0.01;
lambda = 1.0;

logdir=[model 'log/mu=' num2str(mu) '_l=' num2str(lambda)];
load([logdir '/tsp.mat']);

[H, c, fH] = visualize2(points,  orders, [], 'after TSP optimization', logdir, false, true);

S = {}; 
set(H, 'ButtonDownFcn', {@selected2, H});
waitfor(fH);

T  = [];
for it = 1: length(S)
    if length(S{it}) > 0
        T = [T; S{it}, it];
    end
end
T
save([logdir '/tsp.mat'], 'edges', 'orders', 'points', 'D1', 'D2', 'D3', 'D4', 'T')