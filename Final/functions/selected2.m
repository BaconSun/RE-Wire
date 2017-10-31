function selected(ObjectH, EventData, H, clusters, T)
    
    N  = clusters(H == ObjectH);
    NH = H(clusters == N);
    set(NH, 'LineWidth', 4.5);
    set(H(clusters ~= N), 'LineWidth', 1.5);
    NT = T(clusters == N);
    set(NT, 'Visible', 'on');
    set(T(clusters ~= N), 'Visible', 'off');
    
end