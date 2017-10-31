function P = missingcheck(TS0, TS, TT, VS)
    figure;
    plot(TT(1, TS0==0), TT(2, TS0==0), '.r'); hold on;
    %plot(TT(1, :), TT(2, :), '.r'); hold on;
    plot(TT(1, (TS>0)&(TS0==0)), TT(2, (TS>0)&(TS0==0)), '.b'); hold on;
    axis equal;
    
    % check the percentage of cover missing
    for jv = 1:length(VS)
        new_cover = [];
        for kv = 1: size(VS{jv}, 1) 
            new_cover = [new_cover VS{jv}(kv, 2): VS{jv}(kv, 3)];
        end
        new_cover = TS0(new_cover);
        P(jv) = sum(new_cover==0)/length(new_cover);  
    end
end