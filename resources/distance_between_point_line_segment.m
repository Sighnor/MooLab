function [lambda, distance] = distance_between_point_line_segment(p0, p1, p2)
    s = p2 - p1;
    den = dot(s, s);
    lambda = dot((p0 - p1), s) / max(den, 0.001);
    if lambda >= 0 && lambda <= 1 && den >= 0.001
        distance = norm(p0 - (p1 + lambda * s));
    else
        distance_p0p1 = norm(p0 - p1);
        distance_p0p2 = norm(p0 - p2);
        if distance_p0p1 <= distance_p0p2
            lambda = 0;
            distance = distance_p0p1;
        else
            lambda = 1;
            distance = distance_p0p2;
        end
    end
end