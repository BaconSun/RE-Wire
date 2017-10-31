function selected(ObjectH, EventData, H, I)
    
    set(H(ObjectH == H), 'LineWidth', 4.5);
    set(H(ObjectH ~= H), 'LineWidth', 1.5);
    set(I{1}(ObjectH == H), 'Visible', 'on');
    set(I{1}(ObjectH ~= H), 'Visible', 'off');
    set(I{2}(ObjectH == H), 'Visible', 'on');
    set(I{2}(ObjectH ~= H), 'Visible', 'off');
    set(I{3}(ObjectH == H), 'Visible', 'on');
    set(I{3}(ObjectH ~= H), 'Visible', 'off');
end