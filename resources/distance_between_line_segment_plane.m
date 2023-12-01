function [lambda, distance] = distance_between_line_segment_plane(p1, p2, A, B, C, D)
    p = [0, 0, 0]';
    if A ~= 0
        p = [-D / A, 0, 0]';
    elseif B ~= 0
        p = [0, -D / B, 0]';
    elseif C ~= 0
        p = [0, 0, -D / C]';
    end
    dotp1pn = dot(p1 - p, [A, B, C]');
    dotp2pn = dot(p2 - p, [A, B, C]');
    distancep1plane = sign(dotp1pn) * distance_between_point_plane(p1, A, B, C, D);
    distancep2plane = sign(dotp2pn) * distance_between_point_plane(p2, A, B, C, D);
   
    if distancep1plane <= distancep2plane
        lambda = 0;
        distance = distancep1plane;
    else
        lambda = 1;
        distance = distancep2plane;
    end
end