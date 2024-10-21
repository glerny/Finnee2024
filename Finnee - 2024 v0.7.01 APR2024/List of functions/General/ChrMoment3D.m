function M = ChrMoment3D(XY)

if nargin == 1
    n = 4;
end

X = XY{2};
Y = XY{3};
XY = XY{1};

X2Y = [];
for ii = 1:size(XY, 2)
    X2Y(ii) = trapz(Y', XY(:, ii));
end
 M(1, 1) = trapz(X', X2Y');
 M(2, 1) = trapz(X', X'.*X2Y)/M(1, 1);
 M(3, 1) = trapz(X', (X' - M(2, 1)).^2.*X2Y)/M(1, 1);
 M(4, 1) = trapz(X', (X' - M(2, 1)).^3.*X2Y)/M(1, 1);

Y2X = [];
for ii = 1:size(XY, 1)
    Y2X(ii) = trapz(X', XY(ii, :));
end
 M(1, 2) = trapz(Y', Y2X');
 M(2, 2) = trapz(Y', Y.*Y2X')/M(1, 2);
 M(3, 2) = trapz(Y', (Y - M(2, 2)).^2.*Y2X')/M(1, 2);
 M(4, 2) = trapz(Y', (Y - M(2, 2)).^3.*Y2X')/M(1, 2);

end
