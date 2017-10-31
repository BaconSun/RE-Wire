function [x, fval, exitflag] = gurobi_qc(f, intcon, A, b, Aeq, beq, lb, ub, quadcon)

n = size(A, 2);
A = sparse(A);
Aeq = sparse(Aeq);

model.obj = f;
model.vtype = repmat('C', n, 1);
model.vtype(intcon) = 'I';

model.A = [A; Aeq];
model.rhs = [b; beq];
model.sense = [repmat('<', size(A,1), 1); repmat('=', size(Aeq,1), 1)];

model.quadcon = quadcon'; % quad constraints

model.lb = lb;
model.ub = ub;


params.outputflag = 1;
result = gurobi(model, params);


if strcmp(result.status, 'OPTIMAL')
    exitflag = 1;
elseif strcmp(result.status, 'INTERRUPTED')
    if isfield(result, 'x')
        exitflag = 2;
    else
        exitflag = 0;
    end
elseif strcmp(result.status, 'INF_OR_UNBD')
    params.dualreductions = 0;
    result = gurobi(model, params);
    if strcmp(result.status, 'INFEASIBLE')
        exitflag = -2;
    elseif strcmp(result.status, 'UNBOUNDED')
        exitflag = -3;
    else
        exitflag = nan;
    end
else
    exitflag = nan;
end


if isfield(result, 'x')
    x = result.x;
else
    x = nan(n,1);
end

if isfield(result, 'objval')
    fval = result.objval;
else
    fval = nan;
end