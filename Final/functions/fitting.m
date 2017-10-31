function [sp fitpoints save_points]= fitting(F, points, direction, tour, logdir, tol, verbose)

% fit points
fitpoints=[];
for i = 1: length(tour)
    it     = tour(i);
    L      = size(points{it}, 2);
    N      = 30;
    off    = 0;
    if direction(i) == 1
        fitpoints = [fitpoints points{it}(:, 1+off:N:end-off)];
    else
        fitpoints = [fitpoints points{it}(:, end-off:-N:1+off)];
    end
end


order=4;
xarc(1)=0;
knot_num = size(fitpoints, 2);

for i=2:knot_num
    xarc(i) = norm(fitpoints(:,i)-fitpoints(:,i-1))+xarc(i-1);
end

f = spap2(knot_num+order,order,xarc,fitpoints);  % B-spline

%FigHandle1 = figure;
%set(FigHandle1, 'Position', [100, 100, 800, 600]);
%plot3(fitpoints(1,:),fitpoints(2,:),fitpoints(3,:),'o'), hold on
%title('B-spline before smoothing')
%fnplt(f), axis equal


%FigHandle2 = figure;
%set(FigHandle2, 'Position', [100, 100, 800, 600]);
[sp, v, p2]= spaps(xarc,fitpoints,tol);
%plot3(fitpoints(1,:),fitpoints(2,:),fitpoints(3,:),'o'), hold on
%title('B-spline after smoothing');
%fnplt(sp), axis equal, hold on
% sn = fnxtr(sp);
% fnplt(sn), axis equal
%savefig([logdir '/tsp_fitting.fig']);

save_points=[];
pnum=1000;
step=(xarc(end)-xarc(1))/pnum;
for ti=0:step:(xarc(end) - step)
    save_points = [save_points;fnval(sp,ti)'];
end

write_off([logdir '/smooth.off'], save_points, int16.empty(0,3));

% save all the control points
fc = fopen([logdir '/control.txt'],'w');
for it = 1: size(sp.coefs, 2)
  fprintf(fc,'<%f, %f, %f>, 1\n', sp.coefs(1,it), sp.coefs(2,it), sp.coefs(3,it));
end
fclose(fc)


if verbose
figure(F);
    [sp, v, p2]= spaps(xarc,fitpoints,tol);
    %plot3(fitpoints(1,:),fitpoints(2,:),fitpoints(3,:),'o'), hold on
    hold on
    fnplt(sp, 'r'),
    plot3(sp.coefs(1,:), sp.coefs(2,:), sp.coefs(3,:),'k'),
    plot3(sp.coefs(1,:), sp.coefs(2,:), sp.coefs(3,:),'x'),
    title('B-spline with control points')
    axis equal
    view(-36, 32)
    hold on;
end
%savefig([logdir '/tsp_fitting_wcp.fig']);
end