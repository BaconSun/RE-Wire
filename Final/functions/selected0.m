function selected0(ObjectH, EventData, H)
    if get(ObjectH, 'LineWidth') == 4
        flag = 0;
        for it = 1: length(H)
            if get(H(it), 'LineWidth') == 7
                flag = 1

                set(H(it), 'UserData', [get(H(it), 'UserData') find(H==ObjectH)]);
                set(ObjectH, 'UserData', [get(ObjectH, 'UserData') it]);
                set(H(it), 'LineWidth', 1.5);
                set(ObjectH, 'LineWidth', 1.5);
 
                S = {};
                for it = 1: length(H)
                    S{it} = get(H(it), 'UserData');
                end
                assignin('caller', 'S',  S);
                break
            end
        end



        if flag == 0
            set(ObjectH, 'LineWidth', 7);
        end

    else
        set(ObjectH, 'LineWidth', 4);
        for it = 1: length(H)
            if (get(H(it), 'LineWidth') ~= 7) && (H(it)~=ObjectH)
                set(H(it), 'LineWidth', 1.5);
            end
        end
    end
end
