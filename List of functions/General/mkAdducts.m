function myAdducts = mkAdducts(MS_ctr, ThresInAdduct, ThresBtwAdduct)
MaxMZdist = 2.1; % mz

myAdducts = {};
MaxArea = max(MS_ctr(:, 1));

while 1

    if isempty(MS_ctr), break; end

    CurAdduct = MS_ctr(MS_ctr(:, 2) == max(MS_ctr(:, 2)), :);
    if CurAdduct(1,2) < ThresBtwAdduct; break, end
    
    MS_ctr(MS_ctr(:, 2) == max(MS_ctr(:, 2)), :) = [];
    while 1
        IdX = [];
        for ii = 1:height(CurAdduct)
            cIdX = find(abs(CurAdduct(ii, 1) - MS_ctr(:, 1)) <= MaxMZdist & ...
                MS_ctr(:, 2) > ThresInAdduct);
            IdX = unique([IdX; cIdX]);
        end

        if isempty(IdX)
            break
        end

        CurAdduct = [CurAdduct; MS_ctr(IdX, :)];
        CurAdduct = sortrows(CurAdduct, -2);
        MS_ctr(IdX, :) = [];

    end
    myAdducts{end + 1} = CurAdduct;
end

end