function fixednT = findnearest(V, T, k)
    curvedist = [0; sqrt(sum((V(1: end-1, :) - V(2:end, :)) .^ 2, 2))];
    curvedist = cumsum(curvedist);

    nTid  = knnsearch(KDTreeSearcher(T,'distance','euclidean'),V,'k',k);
    
    nT    = T(nTid, :);
    nVid  = knnsearch(KDTreeSearcher(V,'distance','euclidean'),nT,'k',1);
    nVid  = reshape(nVid, [], k);
    [m, mid]  = min(abs(curvedist(nVid) - repmat(curvedist(1: size(V, 1)), 1, k)), [], 2);
    
    fixednTid = nTid(sub2ind(size(nTid), [1: size(V, 1)]', mid));
    fixednT   = T(fixednTid, :);

end