function [f, A, b, Aeq, beq, AvailIdx] = build_mTSP(m, nums, E, D, lamda, avail)
    k    = length(nums);
    L = (k+1)*(k+1);
    
    %% problem
    EdgeCost = zeros(k + 1);
    EdgeCost(1: end-1, 1: end-1) = (D + lamda * E);
    maxCost  = max(max(EdgeCost));
    
    % %% corner egdes has great weights!
    % for ik = 1: size(M, 2)
    %     EdgeCost(M(1, ik), M(2, ik)) = 100000 * maxCost;
    %     EdgeCost(M(2, ik), M(1, ik)) = 100000 * maxCost;
    % end

    EdgeCost = reshape(EdgeCost, [], 1);
    [value index]=sort(EdgeCost);

    avail = [avail; ones(k, 1)];
    AvailIdx = find(avail>0);


    f     = [EdgeCost; zeros(k, 1)];

    %% equality constraint -----------------------------------------------------
    M     = zeros(k + 1);
    
    % 1:outgoing edges from a cluster must be one. (per cluster)
    Aeq1   = sparse(k, L + k);
    beq1   = ones(k, 1);
    for i = 1: k
        M1 = M;
        M1(i, setdiff(1: k+1, i)) = 1;
        M1 = reshape(M1, [], 1);
        Aeq1(i, 1: L) = M1';
    end

    % 2: incoming edges to a cluster must be one.  (per cluster)
    Aeq2   = sparse(k, L + k);
    beq2   = ones(k, 1); 
    for i = 1: k
        M2 = M;
        M2(setdiff(1: k+1, i), i) = 1;
        M2 = reshape(M2, [], 1);
        Aeq2(i, 1: L) = M2';
    end

    % 3: outgoing edges form S must be 1 (or m)
    Aeq3   = sparse(1, L + k);
    beq3   = ones(1, 1) * m;
    M3     = M;
    M3(end, 1:end-1) = 1;
    M3     = reshape(M3, [], 1);
    Aeq3(1, 1: L) = M3';

    % 4: incoming edges to S must be 1 (or m)
    Aeq4   = sparse(1, L + k);
    beq4   = ones(1, 1) * m;
    M4     = M;
    M4(1:end-1, end) = 1;
    M4     = reshape(M4, [], 1);
    Aeq4(1, 1: L) = M4';

    Aeq = [Aeq1; Aeq2; Aeq3; Aeq4];
    beq = [beq1; beq2; beq3; beq4];

    disp('equality constraint complete!');
    size(Aeq)

    %% inequality constraint ----------------------------------------------------
    U   = eye(k);

    % ensure we don't get disconnected cycles:
    % 6:   ui + (k-1-m) * x0i - xi0 <= k-m
    A6  = sparse(k, L + k);
    b6  = ones(k, 1) * (k - m);
    for i = 1: k
        Ui  = U(i, :);

        M61 = M;
        M61(end, i) = 1;  % 0 -> p
        M61 = reshape(M61, [], 1);

        M62 = M;
        M62(i, end) = 1;  % p -> 0
        M62 = reshape(M62, [], 1);

        % constraint:
        M6  = (k - 1 - m) * M61 - M62;
        A6(i, 1: L)       = M6;
        A6(i, L + 1: end) = Ui;
    end

    % ensure we don't get disconnected cycles: (no changes)
    % 7:    up + y0p + (2-k) * yp0 >= 2
    %      -up + (k-2) * yp0 - y0p <= -2
    A7  = sparse(k, L + k);
    b7  = ones(k, 1) * (-2);
    for i = 1: k
        Ui  = U(i, :);

        M71 = M;
        M71(end, i) = 1;  % 0 -> p
        M71 = reshape(M71, [], 1);

        M72 = M;
        M72(i, end) = 1;  % p -> 0
        M72 = reshape(M72, [], 1);

        % constraint:
        M7  = - M71; %(k - 2) * M72 
        A7(i, 1: L)       = M7;
        A7(i, L + 1: end) = -Ui;
    end

    % ensure we don't get disconnected cycles:
    % 8:   up - ul + k * ypl + (k - 2) * ylp <= k - 1
    A8  = sparse(k * k - k, L + k);
    b8  = ones(k * k - k, 1) * (k - m);
    s   = 0;
    for i = 1: k
        for j = 1: k
            
            if i == j
                continue
            end

            s   = s + 1;
            Ui  = U(i, :);
            Uj  = U(j, :);

            M81 = M;
            M81(i, j) = 1;
            M81 = reshape(M81, [], 1);

            M82 = M;
            M82(j, i) = 1;
            M82 = reshape(M82, [], 1);

            % constraint
            M8  = (k + 1 - m) * M81 + (k - 1 - m) * M82;
            A8(s, 1: L)       = M8;
            A8(s, L + 1: end) = Ui - Uj;
        end
    end

    A = [A6; A7; A8];
    b = [b6; b7; b8];
    disp('inequality constraint complete!')
    size(A)
end 