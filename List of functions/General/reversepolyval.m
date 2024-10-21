function x = reversepolyval(x, pSmu)

x0 = x;
options = optimset('Display','none');
[x, ~, exitflag] = fminsearch(@objectivefcn1, x0, options);
count = 1;

while exitflag == 0
    [x, ~, exitflag] = fminsearch(@objectivefcn1, x, options);
    count = count +1;

    if count > 10; break; end

end


    function f = objectivefcn1(a)
        ana = a + polyval(pSmu{1}, a, pSmu{2}, pSmu{3});
        f = abs(x-ana);
    end
end