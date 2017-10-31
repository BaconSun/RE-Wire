function [f, A, b, Aeq, beq, quadcon, AvailIdx] = build_mTSP2(nums, E, D, lamda, avail, v)
    k    = length(nums);
    L    = (k+1)*(k+1);

    %% problem
    EdgeCost = zeros(k + 1);
    EdgeCost(1: end-1, 1: end-1) = (D + lamda * E);
    
    EdgeCost = reshape(EdgeCost, [], 1);
    %[value index]=sort(EdgeCost);

    avail    = [avail; ones(k+1, 1)];
    AvailIdx = find(avail>0);


    f     = [EdgeCost; zeros(k, 1); v]; % xij; uk; m3

    %% equality constraint -----------------------------------------------------
    M0    = zeros(k + 1);
    
    % 1:outgoing edges from a cluster must be one. (per cluster)
    Aeq1   = zeros(k, L + k + 1); %xij, uk, m
    beq1   = ones(k, 1);
    for i = 1: k
        M1 = M0;
        M1(i, setdiff(1: k+1, i)) = 1;
        M1 = reshape(M1, [], 1);
        Aeq1(i, 1: L) = M1';
    end

    % 2: incoming edges to a cluster must be one.  (per cluster)
    Aeq2   = zeros(k, L + k + 1);
    beq2   = ones(k, 1); 
    for i = 1: k
        M2 = M0;
        M2(setdiff(1: k+1, i), i) = 1;
        M2 = reshape(M2, [], 1);
        Aeq2(i, 1: L) = M2';
    end

    % 3: outgoing edges form S must be 1 (or m)
    Aeq3   = zeros(1, L + k + 1);
    beq3   = zeros(1, 1);
    M3     = M0;
    M3(end, 1:end-1) = 1;
    M3     = reshape(M3, [], 1);
    Aeq3(1, 1: L) = M3';
    Aeq3(1, end) = -1;

    % 4: incoming edges to S must be 1 (or m)
    Aeq4   = zeros(1, L + k + 1);
    beq4   = zeros(1, 1);
    M4     = M0;
    M4(1:end-1, end) = 1;
    M4     = reshape(M4, [], 1);
    Aeq4(1, 1: L) = M4';
    Aeq4(1, end) = -1;
    
    Aeq = [Aeq1; Aeq2; Aeq3; Aeq4];
    beq = [beq1; beq2; beq3; beq4];
    disp('equality constraint complete!');

    %% ---- completely different --- %%
    %% inequality
    U = eye(k);
    quadcon = [];
    
    % ensure we don't get disconnected cycles:
    % 6:   ui + (k-1-m) * x0i - xi0 <= k-m (QC)
    % ==>  -m * x0i + ui + (k-1) *x0i - xi0 +m <= k
    Mask = zeros(1, L + k + 1);
    for i = 1: k
        Mui = Mask; Mui(1, L+i) = 1;
        Mm  = Mask; Mm(1,  end) = 1;
        M0i = Mask; M0i(1, getid(k+1, i, k)) = 1;
        Mi0 = Mask; Mi0(1, getid(i, k+1, k)) = 1;
        
        % constraint:
        Qc  = sparse(L+k+1, getid(k+1, i, k), -1, L+k+1, L+k+1);
        q   = Mui + (k-1) * M0i - Mi0 + Mm;
        rhs = k;
        
        quad.Qc  = Qc;
        quad.q   = q';
        quad.rhs = rhs;
        quadcon  = [quadcon; quad];
    end

    % ensure we don't get disconnected cycles: (no changes)
    % 7:    ui + x0i >= 2
    %      -ui - x0i <= -2
    A7  = zeros(k, L + k+1);
    b7  = ones(k, 1) * (-2);
    for i = 1: k
        Ui  = U(i, :);

        M71 = M0;
        M71(end, i) = 1;  % 0 -> p
        M71 = reshape(M71, [], 1);

        % constraint:
        M7  = -M71; %(k - 2) * M72 
        A7(i, 1: L)       = M7;
        A7(i, L + 1: end-1) = -Ui;
    end

    % ensure we don't get disconnected cycles:
    % 8:   ui - uj + (k-m+1) * xij + (k-1-m) * xji <= k - m
    %=>    -m*(xij+xji)+[ui-uj+(k+1)xij+(k-1)xji+m]<=k
    s   = 0;
    for i = 1: k
        for j = 1: k
            if i == j
                continue
            end

            Mui = Mask; Mui(1, L+i) = 1;
            Muj = Mask; Muj(1, L+j) = 1;
            Mm  = Mask; Mm(1,  end) = 1;
            Mji = Mask; Mji(1, getid(j, i, k)) = 1;
            Mij = Mask; Mij(1, getid(i, j, k)) = 1;

            % constraint:
            Qc  = sparse([L+k+1; L+k+1], [getid(j, i, k); getid(i, j, k)], ...
                         [-1, -1], L+k+1, L+k+1);
            q   = Mui - Muj + (k+1) * Mij + (k-1) * Mji + Mm;
            rhs = k;
        
            quad.Qc  = Qc;
            quad.q   = q';
            quad.rhs = rhs;
            quadcon  = [quadcon; quad];
        end
    end

    A = A7;
    b = b7;
    disp('inequality constraint complete!')
    size(A)
    size(quadcon)
end 

function id = getid(i, j, k)
    id = (j - 1) * (k + 1) + i; 
end