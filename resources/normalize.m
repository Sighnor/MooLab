function [output] = normalize(input)
    length = max(norm(input), 0.001);
    output = input / length;
end

