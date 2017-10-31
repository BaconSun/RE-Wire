function check_cluster(ObjectH, EventData, H)
    clusers = cell2mat(get(H, 'UserData'));
    targetc = get(ObjectH, 'UserData');

    set(H(clusers == targetc), 'MarkerSize', 13, 'LineWidth', 4);
    set(H(clusers ~= targetc), 'MarkerSize', 5, 'LineWidth', 1.5);

end