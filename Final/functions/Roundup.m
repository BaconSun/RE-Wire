function [X, E, D, O]=Roundup(C0, E0, D0, nums, mu, l)
% initialization
N    = length(nums);
C = C0; E = E0; D = D0;
X    = -ones(sum(nums), 1);
Xc   = -ones(sum(nums), 1);
ID   = 1: sum(nums);
O    = zeros(sum(nums), 1);

ik   = 0;
while (sum(nums) > N)
    
    lb  = zeros(size(C, 1));
    ub  = ones(size(C, 1));
    A   = [];
    b   = [];
    
    cun = cumsum(nums);
    I   = eye(N);
    
    Aeq = cell2mat(arrayfun(@(x, i) repmat(I(:, i), 1, x), nums, 1:N, 'UniformOutput', false));
    beq = ones(N, 1);
    
    % optimize
    options = optimoptions('quadprog','Display','off');
    [x,fval,exitflag] = quadprog(2*mu*(l*E+D),C,A,b,Aeq,beq,lb,ub,[],options);
    %[x,fval,exitflag] = quadprog(2*mu*E,C,A,b,Aeq,beq,lb,ub,[],options);


    % find closest to 0 or 1 (maybe multiple)
    temp = min(x, 1-x);
    temp(Xc == 1)  = 100; % avoid selected 1 to be chosen again
    [minx, minidx] = min(temp);
    class          = length(find(cun < minidx)) + 1;
    
    if x(minidx) == minx   % the closest is 0, clean this
        X(ID(minidx)) = 0;
        
        nums(class) = nums(class) - 1;
        
        % clean C and E
        C(minidx, :)  = [];
        E(minidx, :)  = [];
        E(:, minidx)  = [];
        D(minidx, :)  = [];
        D(:, minidx)  = [];
        ID(minidx)    = [];
        Xc(minidx)    = [];
        
        % only one value left
        if nums(class) == 1
            rd        = cun(class) - 1;  % last one has been deleted
            
            ik        = ik + 1;
            %sum(nums)
            X(ID(rd)) = 1;
            O(ID(rd)) = ik;
            Xc(rd)    = 1;               % used to avoid multiple access
        end
    else                   % the closest is 1, clean others
        
        if class == 1
            rangec = 1: cun(class);
            rangec(minidx) = [];
        else
            rangec = cun(class - 1) + 1: cun(class);
            rangec(minidx - cun(class - 1)) = [];
        end
        nums(class)= 1;
        
        % clean C, E
        ik            = ik + 1;
        %sum(nums)
        X(ID(rangec)) = 0;
        X(ID(minidx)) = 1;
        O(ID(minidx)) = ik;

        C(rangec, :)  = [];
        E(rangec, :)  = [];
        E(:, rangec)  = [];
        D(rangec, :)  = [];
        D(:, rangec)  = [];
        ID(rangec)    = [];
        Xc(minidx)    = 1;
        Xc(rangec)    = [];
    end
end
% clean un-selected 1s
size(find(X==1))
X(X==-1) = 1;
end