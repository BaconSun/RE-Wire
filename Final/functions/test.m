function test;

    x = 1: 100;
    x = x / 100 * pi;
    y = cos(x);
    z = sin(x);

    h(1) = plot(x, y, 'r-'); hold on;
    h(2) = plot(x, z, 'b-'); hold on;
    xlim([0, 4]);
    ylim([-1, 1]);

    set(h, 'ButtonDownFcn', {@LineSelected, h})

    %w = waitforbuttonpress;
    %if w == 1
    %    set(h1, 'Visible', 'off')
    %end
end



function LineSelected(ObjectH, EventData, H)
    if get(ObjectH, 'LineWidth') == 2.5
        set(ObjectH, 'Visible', 'off');
    else
        set(ObjectH, 'LineWidth', 2.5);
        set(H(H ~= ObjectH), 'LineWidth', 0.5);
    end
end