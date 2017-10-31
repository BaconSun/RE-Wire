clc;clear;tic;

[o, tasks] = setup('');

% modelname='cutecat_syn_n2'
addpath('functions/');
addpath('functions/bspline/');
baseline_path = './line+tsp/';
fid=fopen('result.txt','w+');
% fidn=fopen('result_normalize.txt','w+');
% fid2=fopen('result2.txt','w+');
for i = 1:length(tasks)
    modelname = tasks(i).model
    task   = setup(modelname);
    model  = task.model;
    suffix = task.suffix;
    tol    = task.tol;
    ctol   = task.ctol;
    mu     = task.select(1);
    l      = task.select(2);

    modelpath = ['./tasks/' model '/'];
    camerapath = ['../Models/' modelname '/'];
    logdir    = [modelpath 'log/mu=' num2str(mu) '_l=' num2str(l)];

    % read off file.
    try
        points = read_off([logdir '/final_' model(1:end-1) '.off']);
    catch
        disp(['error to read ' model])
        continue
    end

    vmin = min(points, [], 2);
    vmax = max(points, [], 2);
    vc   = mean(points, 2) % center
    vv   = vmax - vmin
    vlen = norm(vv)

    % down-sample
    disp('-------- down sample -------')
    points = points(:, 1: 100: size(points,2));

    % read image information
    T1 = read_curvepoints([camerapath 'curvepoints1.txt']);
    T2 = read_curvepoints([camerapath 'curvepoints2.txt']);
    T3 = read_curvepoints([camerapath 'curvepoints3.txt']);

    da1 = norm(max(T1)-min(T1));
    da2 = norm(max(T2)-min(T2));
    da3 = norm(max(T3)-min(T3));

    % camera parameters
    try
        M1 = read_camera([camerapath 'camera1.txt']);
        M2 = read_camera([camerapath 'camera2.txt']);
        M3 = read_camera([camerapath 'camera3.txt']);
    catch
        M1 = read_camera2([camerapath 'camera1.txt']);
        M2 = read_camera2([camerapath 'camera2.txt']);
        M3 = read_camera2([camerapath 'camera3.txt']);
    end

    % projection
    disp('-------- projection  -------')
    V1 = project(M1, points);
    V2 = project(M2, points);
    V3 = project(M3, points);

    % find the nearest points
    disp('-------- findnearest -------')
    nT1 = findnearest(V1, T1, 100);
    nT2 = findnearest(V2, T2, 100);
    nT3 = findnearest(V3, T3, 100);

    % distance
    disp('--------  distance   -------')
    dis1 = sum((V1 - nT1).^2, 2);
    dis2 = sum((V2 - nT2).^2, 2);
    dis3 = sum((V3 - nT3).^2, 2);

    dd1 = sqrt(mean(dis1));
    dd2 = sqrt(mean(dis2));
    dd3 = sqrt(mean(dis3));
    dm1 = sqrt(max(dis1));
    dm2 = sqrt(max(dis2));
    dm3 = sqrt(max(dis3));
    dd  = [dd1, dd2, dd3, dm1, dm2, dm3];
    toc;

    % Average on 3 views:     2.244605    4.898622    3.573514    
    % Max on 3 views:     12.727922   16.763055   18.357560
    % Normalized average: 0.0028      0.0038      0.0036
    % Normalized max: 0.0161      0.0129      0.0185
    % Diagonal length:    789.0412    1.2955e+03  991.5891

    % result: values
    fprintf(fid, '%s\n', modelname);
    fprintf(fid, 'Average on 3 views:\t%f\t%f\t%f\n', dd(1), dd(2), dd(3));
    fprintf(fid, 'Maximum on 3 views:\t%f\t%f\t%f\n',     dd(4), dd(5), dd(6));
    fprintf(fid, 'Normalized average:\t%f\t%f\t%f\n', dd(1)/da1, dd(2)/da2, dd(3)/da3);
    fprintf(fid, 'Normalized maximum:\t%f\t%f\t%f\n',     dd(4)/da1, dd(5)/da2, dd(6)/da3);
    fprintf(fid, 'Diagonal length: \t%f\t%f\t%f\n\n',   da1, da2, da3);
    disp('done.');

end

fclose(fid)
% fclose(fidn)
% fclose(fid2)
