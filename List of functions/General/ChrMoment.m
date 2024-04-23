function M = ChrMoment(XY, n)

if nargin == 1
    n = 4;
end

X = XY(:,1);
Y = XY(:,2);

switch n
    case 1
        try
            M(1) = trapz(X, Y);

        catch
            M(1) = NaN;
        end

    case 2
        try
            M(1) = trapz(X, Y);
            M(2) = trapz(X, X.*Y)/M(1);

        catch
            M(1) = NaN;
            M(2) = NaN;
        end

    case 3
        try
            M(1) = trapz(X, Y);
            M(2) = trapz(X, X.*Y)/M(1);
            M(3) = trapz(X, (X - M(2)).^2.*Y)/M(1);

        catch
            M(1) = NaN;
            M(2) = NaN;
            M(3) = NaN;
        end

    case 4
        try
            M(1) = trapz(X, Y);
            M(2) = trapz(X, X.*Y)/M(1);
            M(3) = trapz(X, (X - M(2)).^2.*Y)/M(1);
            M(4) = trapz(X, (X - M(2)).^3.*Y)/M(1);

        catch
            M(1) = NaN;
            M(2) = NaN;
            M(3) = NaN;
            M(4) = NaN;
        end
end
