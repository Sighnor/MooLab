function [lambda1, lambda2, distance] = distance_between_line_segments(p1, p2, q1, q2)
    s1 = p2 - p1;
    s2 = q2 - q1;
    den = (dot(s1, s1) * dot(s2, s2) - dot(s1, s2) * dot(s1, s2));
    lambda1 = (dot(s1, s2) * dot((p1 - q1), s2) - dot(s2, s2) * dot(p1 - q1, s1)) / max(den, 0.001);
    lambda2 = -(dot(s1, s2) * dot((p1 - q1), s1) - dot(s1, s1) * dot(p1 - q1, s2)) / max(den, 0.001);
    if lambda1 >= 0 && lambda1 <= 1 && lambda2 >= 0 && lambda2 <= 1 && den >= 0.001
        distance = norm(p1 + lambda1 * s1 - q1 - lambda2 * s2);
    else
        [lambda_p1l2, distance_p1l2] = distance_between_point_line_segment(p1, q1, q2);
        [lambda_p2l2, distance_p2l2] = distance_between_point_line_segment(p2, q1, q2);
        [lambda_q1l1, distance_q1l1] = distance_between_point_line_segment(q1, p1, p2);
        [lambda_q2l1, distance_q2l1] = distance_between_point_line_segment(q2, p1, p2);
        if distance_p1l2 <= distance_p2l2 && distance_p1l2 <= distance_q1l1 && distance_p1l2 <= distance_q2l1
            lambda1 = 0;
            lambda2 = lambda_p1l2;
            distance = distance_p1l2;
        elseif distance_p2l2 <= distance_q1l1 && distance_p2l2 <= distance_q2l1
            lambda1 = 1;
            lambda2 = lambda_p2l2;
            distance = distance_p2l2;
        elseif distance_q1l1 <= distance_q2l1
            lambda1 = lambda_q1l1;
            lambda2 = 0;
            distance = distance_q1l1;
        else
            lambda1 = lambda_q2l1;
            lambda2 = 1;
            distance = distance_q2l1;
        end
    end
end