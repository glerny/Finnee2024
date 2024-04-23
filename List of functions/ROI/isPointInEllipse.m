function result = isPointInEllipse(x, y, h, k, a, b)
    % Check if a point (x, y) is within an ellipse defined by its center (h, k),
    % semi-major axis length (a), and semi-minor axis length (b).
    
    % Calculate the left-hand side of the ellipse equation
    lhs = ((x - h)^2 / a^2) + ((y - k)^2 / b^2);
    
    % Check if the point is inside or on the boundary of the ellipse
    if lhs <= 1
        result = true;
    else
        result = false;
    end
end