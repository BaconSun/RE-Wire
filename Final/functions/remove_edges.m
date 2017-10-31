function avail = remove_edges(clusters, combined, cc, L)

    disp('remove edges')
    Amat     = ones((L + 1) , (L + 1));
    disp(['before: ' num2str(sum(sum(Amat)))]);

    for ia = 1: (length(combined)-1)
        availmat = ones((L + 1) , (L + 1));
        nstart = combined(ia) + 1;
        nend   = combined(ia+1);

        cstart = clusters(nstart);
        cend   = clusters(nend);

        %availmat(:, cc(cstart)+1: cc(cstart+1))=1;
        %availmat(cc(cend)+1: cc(cend+1), :)=1;

        for c = cstart: (cend-1)
            n_from = cc(c) + 1: cc(c+1);
            n_to   = cc(c+1) +1: cc(c+2);

            if c > cstart
                availmat(n_from, :)    = 0;
                availmat(:, n_from)    = 0;
            end
            if c < (cend-1)
                availmat(n_to, :)      = 0;
                availmat(:, n_to)      = 0;
            end
        end

        for c = cstart: (cend-1)
            n_from = cc(c) + 1: cc(c+1);
            n_to   = cc(c+1) +1: cc(c+2);
            availmat(n_from, n_to) = 1;
            availmat(n_to, n_from) = 1;
        end
        Amat = Amat .* availmat;
    end
    avail = [reshape(Amat, [], 1)];
    disp(['after: ' num2str(sum(sum(Amat)))]);
    %figure; imshow(Amat);


end