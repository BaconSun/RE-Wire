function [X, E, D]=gurobi_intqp(C0, E0, D0, nums, mu, l)
    
    Q   = 2 * mu * (E0 * l + D0);
    N   = length(nums);
    I   = eye(N);
    Aeq = cell2mat(arrayfun(@(x, i) repmat(I(:, i), 1, x), nums, 1:N, 'UniformOutput', false));
    beq = ones(N, 1);

    clear model;
    model.Q = sparse(Q);
    model.A = sparse(Aeq);
    model.obj = C0;
    model.rhs = beq;
    model.sense='=';
    model.vtype  = 'B';
    gurobi_write(model, 'qp.lp')
    r = gurobi(model);
    
    X = r.x;
    E = E0;
    E(X==0, :) = [];
    E(:, X==0) = [];

    D = D0;
    D(X==0, :) = [];
    D(:, X==0) = [];

end