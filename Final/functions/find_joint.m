clc; clear;

model = 'group532/'
source_path = ['../model/' model];
F = importdata([model 'CurveFileName.txt']);
addpath('./joints/');


% load ?
mu=1/4000000; l=100000;
logdir=[model 'slog/mu=' num2str(mu) '_l=' num2str(l)];
load([logdir '/data.mat'])
F = F(X==1);

% V1 curve points
V1     = imread([model 'view1_bd.png']);
V1Edge = segment(V1);
[J1, J1E, J1C] = get_joint(V1Edge);

% V2 curve points
V2     = imread([model 'view2_bd.png']);
V2Edge = segment(V2);
[J2, J2E, J2C] = get_joint(V2Edge);


% load .off points
V1Edge3D = cell(size(V1Edge, 2), 1);
V1EdgeCs = cell(size(V1Edge, 2), 1);
for it = 1 : length(F)
    name   = F{it};
    e      = str2num(name(11:12)) + 1;
    points = read_off([source_path name]);
    points = [points; ones(1, size(points, 2))];
    V1Edge3D{e} = [V1Edge3D{e} points];
    V1EdgeCs{e} = [V1EdgeCs{e} it];
end


% load camera matrix and make pre-projection
matrix = importdata([model 'camera1.txt']);
P1     = matrix(1:3, :);
V1PJ   = cellfun(@(x)project(P1, x), V1Edge3D, 'UniformOutput', false);

matrix = importdata([model 'camera2.txt']);
P2     = matrix(1:3, :);
V2PJ   = cellfun(@(x)project(P2, x), V1Edge3D, 'UniformOutput', false);


% checking joints one by one:
K      = [];
G      = [];
for k  = 1: size(J1, 1)
    if J1C(k) > 1
        dis   = []; idxs = [];
        dis2  = [];
        for j = 1: length(J1E{k})
            ej = J1E{k}(j);
            if ~isempty(V1PJ{ej})
                [v, idx] = min(sqrt((V1PJ{ej}(:, 1) - J1(k, 2)).^2 ...
                                  + (V1PJ{ej}(:, 2) - J1(k, 1)).^2));  % Note that x/y reverse.
                dis = [dis v]; idxs = [idxs idx];
            
                v2  = sqrt((V2PJ{ej}(idx, 1) - J2(:, 2)).^2 ...
                         + (V2PJ{ej}(idx, 2) - J2(:, 1)).^2);          % Note that x/y reverse.
                dis2 = [dis2  v2];
            end
        end

        D1     = min(dis);
        dis2   = min(dis2, [], 2);
        [D2 g] = min(dis2);
        
        if D2 / (D1 + 1e-7) > 20.0  % Threshold, you can change it.
            disp('fake');
        else
            disp('real');
            K = [K k];
            G = [G g];
        end

    end
end

% visualization
figure(20);
for i = 1: size(V1Edge, 2)
    plot(V1Edge{1, i}(:, 2), -V1Edge{1, i}(:, 1), 'r.'); hold on;
end
for r = 1: length(K)
    plot(J1(K(r), 2), -J1(K(r), 1), 'bx', 'MarkerSize', 40); hold on;
end
axis equal;
%saveas(gcf,'real_joint.png')

figure(30);
for i = 1: size(V2Edge, 2)
    plot(V2Edge{1, i}(:, 2), -V2Edge{1, i}(:, 1), 'r.'); hold on;
end
for r = 1: length(G)
    plot(J2(G(r), 2), -J2(G(r), 1), 'bx', 'MarkerSize', 40); hold on;
end
axis equal;
%saveas(gcf,'real_joint2.png')

save([logdir '/joints.mat'], 'V1EdgeCs', 'J1E', 'K');
disp('done.')



