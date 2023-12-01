function [distance] = distance_between_point_plane(p, A, B, C, D)
    distance = abs(A * p(1) + B * p(2) + C * p(3) + D) / sqrt(A * A + B * B + C * C);
end