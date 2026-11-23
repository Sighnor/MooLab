function [output] = quat_to_R(input)
    num = length(input);
    output = zeros(3, 3 * num);
    for i = 1:num
        theta = 180 / pi * (2 * acos(input(1, i)));
        omega = normalize(input(2:4, i));
        output(:, 3 * i - 2:3 * i) = Rodrigues(theta, omega);
    end
end

