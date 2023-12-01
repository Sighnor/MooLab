function [output] = Rodrigues(theta,omega)
    my_theta = theta * pi / 180;
    v = 1 - cos(my_theta);
    c = cos(my_theta);
    s = sin(my_theta);  
    w1 = omega(1);
    w2 = omega(2);
    w3 = omega(3);  
    
    output = [w1^2 * v + c,w1 * w2 * v - w3 * s,w1 * w3 * v + w2 * s;
                 w1 * w2 * v + w3 * s,w2^2 * v + c,w2 * w3 * v - w1 * s;
                 w1 * w3 * v - w2 * s,w2 * w3 * v + w1 * s,w3^2 * v + c];
end

