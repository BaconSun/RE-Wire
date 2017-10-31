function selected4(ObjectH, EventData, H, VD)
    N  = find(ObjectH==H);
    NA = VD{N};
    set(H(NA), 'LineWidth', 4.5);
    set(H(setdiff(1: length(H), NA)), 'LineWidth', 1.5);
   
end