function [Joints, JointEdges, JointCount] = get_joint(nedgelist)
    % return joint points from edge lists
    
    EndPoints=zeros(2 * size(nedgelist, 2), 2);
    for i=1:size(nedgelist, 2)
        Vertices=nedgelist{1,i};
        EndPoints(2 * i - 1, :) = Vertices(1,:);
        EndPoints(2 * i, :)     = Vertices(end,:);
    end
    Edges = reshape(repmat((1: size(nedgelist, 2))', 1, 2).', 1, []);

    % find joints directly using unique
    [Joints, ri, ui] = unique(EndPoints, 'rows');
    JointEdges = cell(length(Joints), 1);
    for k = 1: length(ui)
        JointEdges{ui(k)} = [JointEdges{ui(k)}; Edges(k)];
    end
    JointCount  = cellfun(@length, JointEdges);
end

