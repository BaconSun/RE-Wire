function selected(ObjectH, EventData, H)
    if get(ObjectH, 'LineWidth') == 4
        set(ObjectH, 'Visible', 'off');
        for it = 1: length(H)
            X(it) = strcmp(get(H(it), 'Visible'), 'on');
        end

        assignin('caller', 'X',  X);
    else
        set(ObjectH, 'LineWidth', 4);
        set(H(H ~= ObjectH), 'LineWidth', 1.5);
    end
end