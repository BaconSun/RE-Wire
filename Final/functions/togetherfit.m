function [save_points dd]=togetherfit(model, sp, tol, ctol, pnum, iters, verbose)
    
    showIter = verbose;

    % final fitting with images
    % read image points!
    tic
    T1 = read_curvepoints([model 'curvepoints1.txt']);
    T2 = read_curvepoints([model 'curvepoints2.txt']);
    T3 = read_curvepoints([model 'curvepoints3.txt']);

    % camera parameters
    try
        M1 = read_camera([model 'camera1.txt']);
        M2 = read_camera([model 'camera2.txt']);
        M3 = read_camera([model 'camera3.txt']);
    catch
        M1 = read_camera2([model 'camera1.txt']);
        M2 = read_camera2([model 'camera2.txt']);
        M3 = read_camera2([model 'camera3.txt']);
    end
    % points on B-spline
    N = round(pnum / 100);     % # of samples


    % inital control points
    P  = sp.coefs;            % 3 x K  initial control points
    K  = size(P, 2);
    P2 = reshape(P, [], 1);   % 3K x 1  

    % B-spline matrix
    C  = bspline_basismatrix(sp.order, sp.knots, ...
         linspace(sp.knots(sp.order), sp.knots(end-sp.order+1), N)); % N x K
    C2 = zeros(3 * N, 3 * K); % 3N x 3K
    C2(1:3:end, 1:3:end) = C; C2(2:3:end, 2:3:end) = C; C2(3:3:end, 3:3:end) = C;

    % smooth factor matrix
    G  = eye(N-2, N);
    G  = G - 2 * [G(:, end) G(:, 1:end-1)] + [G(:, end-1: end) G(:, 1:end-2)];
    G2 = zeros(3 * (N-2), 3 * N);
    G2(1:3:end, 1:3:end) = G; G2(2:3:end, 2:3:end) = G; G2(3:3:end, 3:3:end) = G;
    S2 = G2 * C2;    % 3(N-2) x 3N

    % start iteration!!!
    for it = 1: (iters+1)
        if iters == 0
            break
        end
        
        
        points = reshape(C2 * P2, 3, [])';  % 3D curve points
        
        % projection
        V1 = project(M1, points);
        V2 = project(M2, points);
        V3 = project(M3, points);

        % find the nearest points
%         nT1 = findnearest(V1, T1, 50);
%         nT2 = findnearest(V2, T2, 50);
%         nT3 = findnearest(V3, T3, 50);
        
        nT1 = findnearest(V1, T1, 25);
        nT2 = findnearest(V2, T2, 25);
        nT3 = findnearest(V3, T3, 25);
        
        % distance
        dis1 = sum((V1 - nT1).^2, 2);
        dis2 = sum((V2 - nT2).^2, 2);
        dis3 = sum((V3 - nT3).^2, 2);
        dd1 = sqrt(mean(dis1));
        dd2 = sqrt(mean(dis2));
        dd3 = sqrt(mean(dis3));
        dm1 = sqrt(max(dis1));
        dm2 = sqrt(max(dis2));
        dm3 = sqrt(max(dis3));
        dd  = [dd1, dd2, dd3, dm1, dm2, dm3]
        
        %return
        %nT1  = T1(knnsearch(KDTreeSearcher(T1,'distance','euclidean'),V1,'k',1), :);
        %nT2  = T2(knnsearch(KDTreeSearcher(T2,'distance','euclidean'),V2,'k',1), :);
        %nT3  = T3(knnsearch(KDTreeSearcher(T3,'distance','euclidean'),V3,'k',1), :);

        if showIter
            if (it == (iters+1) | (it == 1))
                fh(it) = figure;
                set(fh(it), 'Position', [0 0 1000 1000]);
               
                subplot(221);
                %title(['results@Iter=' num2str(it)])
                plot3(points(:, 1), points(:, 2), -points(:, 3), '.m'); axis equal;
                
                subplot(222);
                plot(T1(:, 1),  -T1(:, 2),  '.r'); hold on;
                plot(nT1(:, 1), -nT1(:, 2), '.g'); hold on;
                plot(V1(:, 1),  -V1(:, 2),  '.b'); hold on;
                axis equal;legend('image', 'nearest', 'spline');

                subplot(223);
                plot(T2(:, 1),  -T2(:, 2),  '.r'); hold on;
                plot(nT2(:, 1), -nT2(:, 2), '.g'); hold on;
                plot(V2(:, 1),  -V2(:, 2),  '.b'); hold on;
                axis equal;legend('image', 'nearest', 'spline');

                subplot(224);
                plot(T3(:, 1),  -T3(:, 2),  '.r'); hold on;
                plot(nT3(:, 1), -nT3(:, 2), '.g'); hold on;
                plot(V3(:, 1),  -V3(:, 2),  '.b'); hold on;
                axis equal;legend('image', 'nearest', 'spline');

                text(0.5, 1,['\bf results@Iter=' num2str(it-1)], ...
                    'HorizontalAlignment','center', ...
                    'VerticalAlignment', 'top', ...
                    'FontSize',14)
            end
        end

        if it > iters
            break
        end

        % build optimization problem
        A1 = repmat(M1(1:2, 1:3), size(nT1, 1), 1) - reshape(nT1', [], 1) * M1(3, 1:3);  % (2xN) x 3
        A2 = repmat(M2(1:2, 1:3), size(nT2, 1), 1) - reshape(nT2', [], 1) * M2(3, 1:3);  % (2xN) x 3
        A3 = repmat(M3(1:2, 1:3), size(nT3, 1), 1) - reshape(nT3', [], 1) * M3(3, 1:3);  % (2xN) x 3
        A  = [A1; A2; A3];

        % smoothness fatcors
        b1 = reshape(nT1', [], 1) * M1(3, 4) - repmat(M1(1:2, 4), size(nT1, 1), 1);       % (2xN) x 1
        b2 = reshape(nT2', [], 1) * M2(3, 4) - repmat(M2(1:2, 4), size(nT2, 1), 1);       % (2xN) x 1
        b3 = reshape(nT3', [], 1) * M3(3, 4) - repmat(M3(1:2, 4), size(nT3, 1), 1);       % (2xN) x 1
        b  = [b1; b2; b3];

        C3 = [C2; C2; C2];

        B  = zeros(6 * N, 3 * K);
        for j = 1: (3*N)
            B(2*j-1: 2*j, :) = A(2*j-1: 2*j, :) * C3(3*j-2: 3*j, :);
        end

        BS = [B; ctol * S2];
        bs = [b; zeros(3*(N-2), 1)];

        % solve linear function
        P2 = BS \ bs;  % <---------- new control points

    end
    elapsedTime3 = toc

    % save the final curve
    % delete([logdir '/Fitting_info.txt']);
    % fileID = fopen([logdir '/Fitting_info.txt'],'w');
    % fprintf(fileID,'model\t %s\n',model);
    % fprintf(fileID,'tol=1/%f, ctol=%f\n',1/tol,ctol);
    % fprintf(fileID,'Fitting 3D PointsTime:\t %f sec\n' ,elapsedTime2);
    % fprintf(fileID,'Fitting Views Time:\t %f sec\n',elapsedTime3);
    % fclose(fileID);
    % disp('Save Fitting information');

    C  = bspline_basismatrix(sp.order, sp.knots, ...
         linspace(sp.knots(sp.order), sp.knots(end-sp.order+1), pnum)); % N x K
    P  = reshape(P2, 3, []);
    save_points = P * C';
    %write_off([logdir '/final_' model(1:end-1) '.off'], save_points, int16.empty(0,3));
    %disp('All done.');

end
