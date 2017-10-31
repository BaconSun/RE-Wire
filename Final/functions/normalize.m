function D = normalize(Dir)
    if norm(Dir) > 0
        D = Dir / norm(Dir);
    else
        D = Dir;
    end
end