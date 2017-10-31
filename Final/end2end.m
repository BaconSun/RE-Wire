clc; clear;

% The model's name
model = 'human';

warning(['OFF']);
addpath('./functions')
addpath('functions/bspline/');
addpath('./matlabExe/')
if ~exist('results')
    mkdir('results');
end

if ~exist(['results/' model]);
    mkdir(['results/' model]);
end

% load setting
[task tasks options]   = setup(model);

options.given_m     = true;

fileID = fopen(['results/' model '/time_info.txt'], 'w');
tic;
% ------------------------------ %
% reconstruction (with exe)
RecoCurveCore('matlabExe/', '../Models/', 'tasks/',options.name, options)
rec_time = toc; tic

% preprocess
preprocess(task, false);
pre_time = toc; tic;

% select curves + QO
filter_curves(task, options, false);
fil_time = toc; tic;

%pick_corner(task);

% TSP optimization
runTSP(task, options.given_m, false);
tsp_time = toc; tic;

% Bspline+Image Fitting
finalfit(task, false);
fit_time = toc;

fprintf(fileID, 'reconstruction time:\t %f sec\n',   rec_time);
fprintf(fileID, 'preprocessing time:\t %f sec\n',    pre_time);
fprintf(fileID, 'filtering curves time:\t %f sec\n', fil_time);
fprintf(fileID, 'running TSP time:\t %f sec\n',      tsp_time);
fprintf(fileID, 'spline fitting time:\t %f sec\n',   fit_time);
fclose(fileID)

disp('All Done!!');