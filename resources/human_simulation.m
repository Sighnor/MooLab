num_bones = 22;
num_links = 21;
times = 3000;

array_m = zeros(1, num_bones);
for i = 1:num_bones
    array_m(i) = 1;
end

array_I = zeros(3, 3 * num_bones);
for i = 1:num_bones
    array_I(:, 3 * i - 2:3 * i) = eye(3);
end

array_links = [1, 2, 3, 4, 1, 6, 7, 8, 1, 10, 11, 12, 13, 12, 15, 16, 17, 12, 19, 20, 21; 
               2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22];

array_rl = zeros(2, num_bones);
array_rl(1:2, 1) = [0.5 * 0.0738, 0.0738]';
array_rl(1:2, 2) = [0.5 * 0.1, 0.435]';
array_rl(1:2, 3) = [0.5 * 0.1, 0.4237]';
array_rl(1:2, 4) = [0.5 * 0.1, 0.173]';
array_rl(1:2, 5) = [0.5 * 0.05, 0.05]';
array_rl(1:2, 6) = [0.5 * 0.1, 0.435]';
array_rl(1:2, 7) = [0.5 * 0.1, 0.4237]';
array_rl(1:2, 8) = [0.5 * 0.1, 0.173]';
array_rl(1:2, 9) = [0.5 * 0.05, 0.05]';
array_rl(1:2, 10) = [0.5 * 0.12, 0.1259]';
array_rl(1:2, 11) = [0.5 * 0.12, 0.1234]';
array_rl(1:2, 12) = [0.5 * 0.12, 0.2583]';
array_rl(1:2, 13) = [0.5 * 0.10, 0.1177]';
array_rl(1:2, 14) = [0.5 * 0.15, 0.25]';
array_rl(1:2, 15) = [0.5 * 0.1128, 0.2583]';
array_rl(1:2, 16) = [0.5 * 0.1, 0.33]';
array_rl(1:2, 17) = [0.5 * 0.1, 0.252]';
array_rl(1:2, 18) = [0.5 * 0.1, 0.2]';
array_rl(1:2, 19) = [0.5 * 0.1128, 0.2583]';
array_rl(1:2, 20) = [0.5 * 0.1, 0.33]';
array_rl(1:2, 21) = [0.5 * 0.1, 0.252]';
array_rl(1:2, 22) = [0.5 * 0.1, 0.2]';

array_vt = zeros(3, num_bones);
array_wt = zeros(3, num_bones);

array_xt = zeros(3, num_bones);

% array_xt(:, 1) = [0, 0, 0]';
% array_xt(:, 2) = [-1, 0, 0]';
% array_xt(:, 3) = [-2, 0, 0]';
% array_xt(:, 4) = [-3, 0, 0]';
% array_xt(:, 5) = [-4, 0, 0]';
% array_xt(:, 6) = [1, 0, 0]';
% array_xt(:, 7) = [2, 0, 0]';
% array_xt(:, 8) = [3, 0, 0]';
% array_xt(:, 9) = [4, 0, 0]';
% array_xt(:, 10) = [0, 1, 0]';
% array_xt(:, 11) = [0, 2, 0]';
% array_xt(:, 12) = [0, 3, 0]';
% array_xt(:, 13) = [0, 4, 0]';
% array_xt(:, 14) = [0, 5, 0]';
% array_xt(:, 15) = [-1, 3, 0]';
% array_xt(:, 16) = [-2, 3, 0]';
% array_xt(:, 17) = [-3, 3, 0]';
% array_xt(:, 18) = [-4, 3, 0]';
% array_xt(:, 19) = [1, 3, 0]';
% array_xt(:, 20) = [2, 3, 0]';
% array_xt(:, 21) = [3, 3, 0]';
% array_xt(:, 22) = [4, 3, 0]';

array_xt(:, 1) = [0, 0, 0]';
array_xt(:, 2) = [0.5 * 0.211, -0.5 * array_rl(2, 1) - 0.5 * array_rl(2, 2), 0]';
array_xt(:, 3) = [0.5 * 0.211, -0.5 * array_rl(2, 1) - 1 * array_rl(2, 2) - 0.5 * array_rl(2, 3), 0]';
array_xt(:, 4) = [0.5 * 0.211, -0.5 * array_rl(2, 1) - 1 * array_rl(2, 2) - 1 * array_rl(2, 3) - 0.5 * array_rl(2, 4), 0]';
array_xt(:, 5) = [0.5 * 0.211, -0.5 * array_rl(2, 1) - 1 * array_rl(2, 2) - 1 * array_rl(2, 3) - 1 * array_rl(2, 4) - 0.5 * array_rl(2, 5), 0]';
array_xt(:, 6) = [-0.5 * 0.211, -0.5 * array_rl(2, 1) - 0.5 * array_rl(2, 6), 0]';
array_xt(:, 7) = [-0.5 * 0.211, -0.5 * array_rl(2, 1) - 1 * array_rl(2, 6) - 0.5 * array_rl(2, 7), 0]';
array_xt(:, 8) = [-0.5 * 0.211, -0.5 * array_rl(2, 1) - 1 * array_rl(2, 6) - 1 * array_rl(2, 7) - 0.5 * array_rl(2, 8), 0]';
array_xt(:, 9) = [-0.5 * 0.211, -0.5 * array_rl(2, 1) - 1 * array_rl(2, 6) - 1 * array_rl(2, 7) - 1 * array_rl(2, 8) - 0.5 * array_rl(2, 9), 0]';
array_xt(:, 10) = [0, 0.5 * array_rl(2, 1) + 0.5 * array_rl(2, 10), 0]';
array_xt(:, 11) = [0, 0.5 * array_rl(2, 1) + 1 * array_rl(2, 10) + 0.5 * array_rl(2, 11), 0]';
array_xt(:, 12) = [0, 0.5 * array_rl(2, 1) + 1 * array_rl(2, 10) + 1 * array_rl(2, 11) + 0.5 * array_rl(2, 12), 0]';
array_xt(:, 13) = [0, 0.5 * array_rl(2, 1) + 1 * array_rl(2, 10) + 1 * array_rl(2, 11) + 1 * array_rl(2, 12) + 0.5 * array_rl(2, 13), 0]';
array_xt(:, 14) = [0, 0.5 * array_rl(2, 1) + 1 * array_rl(2, 10) + 1 * array_rl(2, 11) + 1 * array_rl(2, 12) + 1 * array_rl(2, 13) + 0.5 * array_rl(2, 14), 0]';
array_xt(:, 15) = [array_rl(1, 12) + array_rl(1, 15), 0.5 * array_rl(2, 1) + 1 * array_rl(2, 10) + 1 * array_rl(2, 11) + 0.5 * array_rl(2, 15), 0]';
array_xt(:, 16) = [array_rl(1, 12) + 1 * 0.1128 + 0.5 * array_rl(2, 16), 0.5 * array_rl(2, 1) + 1 * array_rl(2, 10) + 1 * array_rl(2, 11) + 0.5 * array_rl(2, 15), 0]';
array_xt(:, 17) = [array_rl(1, 12) + 1 * 0.1128 + 1 * array_rl(2, 16) + 0.5 * array_rl(2, 17), 0.5 * array_rl(2, 1) + 1 * array_rl(2, 10) + 1 * array_rl(2, 11) + 0.5 * array_rl(2, 15), 0]';
array_xt(:, 18) = [array_rl(1, 12) + 1 * 0.1128 + 1 * array_rl(2, 16) + 1 * array_rl(2, 17) + 0.5 * array_rl(2, 18), 0.5 * array_rl(2, 1) + 1 * array_rl(2, 10) + 1 * array_rl(2, 11) + 0.5 * array_rl(2, 15), 0]';
array_xt(:, 19) = [-array_rl(1, 12) - array_rl(1, 19), 0.5 * array_rl(2, 1) + 1 * array_rl(2, 10) + 1 * array_rl(2, 11) + 0.5 * array_rl(2, 19), 0]';
array_xt(:, 20) = [-array_rl(1, 12) - 1 * 0.1128 - 0.5 * array_rl(2, 20), 0.5 * array_rl(2, 1) + 1 * array_rl(2, 10) + 1 * array_rl(2, 11) + 0.5 * array_rl(2, 19), 0]';
array_xt(:, 21) = [-array_rl(1, 12) - 1 * 0.1128 - 1 * array_rl(2, 20) - 0.5 * array_rl(2, 21), 0.5 * array_rl(2, 1) + 1 * array_rl(2, 10) + 1 * array_rl(2, 11) + 0.5 * array_rl(2, 19), 0]';
array_xt(:, 22) = [-array_rl(1, 12) - 1 * 0.1128 - 1 * array_rl(2, 20) - 1 * array_rl(2, 21) - 0.5 * array_rl(2, 22), 0.5 * array_rl(2, 1) + 1 * array_rl(2, 10) + 1 * array_rl(2, 11) + 0.5 * array_rl(2, 19), 0]';

% array_xt = textread('X0.txt')';

array_Rt = zeros(3, 3 * num_bones);
for i = 1:num_bones
    array_Rt(:, 3 * i - 2:3 * i) = eye(3);
end

for i = 16:18
    array_Rt(:, 3 * i - 2:3 * i) = Rodrigues(-90, [0, 0, 1]);
end
for i = 20:22
    array_Rt(:, 3 * i - 2:3 * i) = Rodrigues(90, [0, 0, 1]);
end

% array_Rt = quat_to_R(textread('R0.txt')');

array_rt = zeros(3 * 2, num_links);

% for j = 1:num_links
%     array_rt(1:3, j) = 0.5 * (array_xt(:, array_links(2, j)) - array_xt(:, array_links(1, j)));
%     array_rt(4:6, j) = 0.5 * (array_xt(:, array_links(1, j)) - array_xt(:, array_links(2, j)));
% end

array_rt(1:3, 1) = [0.5 * 0.211, -0.5 * array_rl(2, 1), 0];
array_rt(4:6, 1) = [0, 0.5 * array_rl(2, 2), 0];
array_rt(1:3, 2) = [0, -0.5 * array_rl(2, 2), 0];
array_rt(4:6, 2) = [0, 0.5 * array_rl(2, 3), 0];
array_rt(1:3, 3) = [0, -0.5 * array_rl(2, 3), 0];
array_rt(4:6, 3) = [0, 0.5 * array_rl(2, 4), 0];
array_rt(1:3, 4) = [0, -0.5 * array_rl(2, 4), 0];
array_rt(4:6, 4) = [0, 0.5 * array_rl(2, 5), 0];
array_rt(1:3, 5) = [-0.5 * 0.211, -0.5 * array_rl(2, 1), 0];
array_rt(4:6, 5) = [0, 0.5 * array_rl(2, 6), 0];
array_rt(1:3, 6) = [0, -0.5 * array_rl(2, 6), 0];
array_rt(4:6, 6) = [0, 0.5 * array_rl(2, 7), 0];
array_rt(1:3, 7) = [0, -0.5 * array_rl(2, 7), 0];
array_rt(4:6, 7) = [0, 0.5 * array_rl(2, 8), 0];
array_rt(1:3, 8) = [0, -0.5 * array_rl(2, 8), 0];
array_rt(4:6, 8) = [0, 0.5 * array_rl(2, 9), 0];
array_rt(1:3, 9) = [0, 0.5 * array_rl(2, 1), 0];
array_rt(4:6, 9) = [0, -0.5 * array_rl(2, 10), 0];
array_rt(1:3, 10) = [0, 0.5 * array_rl(2, 10), 0];
array_rt(4:6, 10) = [0, -0.5 * array_rl(2, 11), 0];
array_rt(1:3, 11) = [0, 0.5 * array_rl(2, 11), 0];
array_rt(4:6, 11) = [0, -0.5 * array_rl(2, 12), 0];
array_rt(1:3, 12) = [0, 0.5 * array_rl(2, 12), 0];
array_rt(4:6, 12) = [0, -0.5 * array_rl(2, 13), 0];
array_rt(1:3, 13) = [0, 0.5 * array_rl(2, 13), 0];
array_rt(4:6, 13) = [0, -0.5 * array_rl(2, 14), 0];
array_rt(1:3, 14) = [array_rl(1, 12), 0, 0];
array_rt(4:6, 14) = [-array_rl(1, 15), 0, 0];
array_rt(1:3, 15) = [array_rl(1, 15), 0, 0];
array_rt(4:6, 15) = [-0.5 * array_rl(2, 16), 0, 0];
array_rt(1:3, 16) = [0.5 * array_rl(2, 16), 0, 0];
array_rt(4:6, 16) = [-0.5 * array_rl(2, 17), 0, 0];
array_rt(1:3, 17) = [0.5 * array_rl(2, 17), 0, 0];
array_rt(4:6, 17) = [-0.5 * array_rl(2, 18), 0, 0];
array_rt(1:3, 18) = [-array_rl(1, 12), 0, 0];
array_rt(4:6, 18) = [array_rl(1, 19), 0, 0];
array_rt(1:3, 19) = [-array_rl(1, 19), 0, 0];
array_rt(4:6, 19) = [0.5 * array_rl(2, 20), 0, 0];
array_rt(1:3, 20) = [-0.5 * array_rl(2, 20), 0, 0];
array_rt(4:6, 20) = [0.5 * array_rl(2, 21), 0, 0];
array_rt(1:3, 21) = [-0.5 * array_rl(2, 21), 0, 0];
array_rt(4:6, 21) = [0.5 * array_rl(2, 22), 0, 0];

array_delta = zeros(3, num_links);

h = 1 / 60;

x = zeros(3 * times, num_bones);
r = zeros(3 * times, 3 * num_bones);

for i = 1:times
    array_It = zeros(3, 3 * num_bones);
    F = zeros(6 * num_bones, 1);
    vt = zeros(6 * num_bones, 1);
    M = zeros(6 * num_bones, 6 * num_bones);
    J = zeros(3 * num_links, 6 * num_bones);
    delta = zeros(3 * num_links, 1);

    for j = 1:num_bones
        array_It(:, 3 * j - 2:3 * j) = array_Rt(:, 3 * j - 2:3 * j) * array_I(:, 3 * j - 2:3 * j) * array_Rt(:, 3 * j - 2:3 * j)';

        F(6 * j - 5:6 * j - 3) = [0, 9.8 * (-array_m(j)), 0]';
        F(6 * j - 2:6 * j) = [0, 0, 0]' - cross(array_wt(:, j), array_It(:, 3 * j - 2:3 * j) * array_wt(:, j));

        vt(6 * j - 5:6 * j - 3) = array_vt(:, j);
        vt(6 * j - 2:6 * j) = array_wt(:, j);

        M(6 * j - 5:6 * j - 3, 6 * j - 5:6 * j - 3) = 1 / h * array_m(j) * eye(3);
        M(6 * j - 2:6 * j, 6 * j - 2:6 * j) = 1 / h * array_It(:, 3 * j - 2:3 * j);
    end

%     F(6 * 5 - 5:6 * 5 - 3) = [0, 0.33 * 2, 0]';
%     F(6 * 9 - 5:6 * 9 - 3) = [0, 0.33 * 2, 0]';
%     F(6 * 14 - 5:6 * 14 - 3) = F(6 * 14 - 5:6 * 14 - 3) + [0, 9.8 * 22.5, 0]';
%     F(6 * 18 - 5:6 * 18 - 3) = F(6 * 18 - 5:6 * 18 - 3) + [-5, 0, 0]';

    if i < 50
        F(6 * 14 - 5:6 * 14 - 3) = F(6 * 14 - 5:6 * 14 - 3) + [10, 0, 0]';
    end

    for j = 1:num_bones
        eps = 1e-2;
        p1 = array_xt(:, j) + array_Rt(:, 3 * j - 2:3 * j) * [0, 0.5 * array_rl(2, j) - array_rl(1, j), 0]';
        p2 = array_xt(:, j) - array_Rt(:, 3 * j - 2:3 * j) * [0, 0.5 * array_rl(2, j) - array_rl(1, j), 0]';
        for k = j + 1:num_bones
            q1 = array_xt(:, k) + array_Rt(:, 3 * k - 2:3 * k) * [0, 0.5 * array_rl(2, k) - array_rl(1, k), 0]';
            q2 = array_xt(:, k) - array_Rt(:, 3 * k - 2:3 * k) * [0, 0.5 * array_rl(2, k) - array_rl(1, k), 0]';
            [lambda1, lambda2, distance] = distance_between_line_segments(p1, p2, q1, q2);
            if distance < array_rl(1, j) + array_rl(1, k) - eps
                mag_F = 1000 * (array_rl(1, j) + array_rl(1, k) - eps - distance);
                mag_tao = 1000 * (array_rl(1, j) + array_rl(1, k) - eps - distance);
                p = p1 + lambda1 * (p2 - p1);
                q = q1 + lambda2 * (q2 - q1);
                F(6 * j - 5:6 * j - 3) = F(6 * j - 5:6 * j - 3) + mag_F * normalize(p - q);
                F(6 * j - 2:6 * j) = F(6 * j - 2:6 * j) + mag_tao * cross(p - array_xt(:, j), normalize(p - q));
                F(6 * k - 5:6 * k - 3) = F(6 * k - 5:6 * k - 3) + mag_F * normalize(q - p);
                F(6 * k - 2:6 * k) = F(6 * k - 2:6 * k) + mag_tao * cross(q - array_xt(:, k), normalize(q - p));
            end
        end
        [lambda, distance] = distance_between_line_segment_plane(p1, p2, 0, 1, 0, 1.1);
        if distance < array_rl(1, j) - eps
            magF = 1000 * (array_rl(1, j) - eps - distance);
            mag_tao = 1000 * (array_rl(1, j) - eps - distance);
            p = p1 + lambda * (p2 - p1);
            F(6 * j - 5:6 * j - 3) = F(6 * j - 5:6 * j - 3) + mag_F * [0, 1, 0]';
            F(6 * j - 2:6 * j) = F(6 * j - 2:6 * j) + mag_tao * cross(p - array_xt(:, j), [0, 1, 0]');
        end
    end

    total_delta = 0;
    for j = 1:num_links
        id1 = array_links(1, j);
        id2 = array_links(2, j);
        array_delta(:, j) = array_xt(:, id1) + array_rt(1:3, j) - array_xt(:, id2) - array_rt(4:6, j);
        total_delta = total_delta + norm(array_xt(:, id1) + array_rt(1:3, j) - array_xt(:, id2) - array_rt(4:6, j));
    end
    total_delta
    
    for j = 1:num_links
        id1 = array_links(1, j);
        id2 = array_links(2, j);
        J(3 * j - 2:3 * j, 6 * id1 - 5:6 * id1 - 3) = eye(3);
        J(3 * j - 2:3 * j, 6 * id1 - 2:6 * id1) = -cross_mat(array_rt(1:3, j));
        J(3 * j - 2:3 * j, 6 * id2 - 5:6 * id2 - 3) = -eye(3);
        J(3 * j - 2:3 * j, 6 * id2 - 2:6 * id2) = cross_mat(array_rt(4:6, j));
        delta(3 * j - 2:3 * j, 1) = -array_delta(:, j) / h;
    end

    A = [M, - J';
         J, zeros(3 * num_links, 3 * num_links)];
    B = [M * vt + F;delta];

    vt_1 = inv(A) * B;
    
    for j = 1:num_bones
        array_vt(:, j) = vt_1(6 * j - 5:6 * j - 3);
        array_wt(:, j) = vt_1(6 * j - 2:6 * j);
        array_xt(:, j) = array_vt(:, j) * h + array_xt(:, j);
        array_Rt(:, 3 * j - 2:3 * j) = Rodrigues(180 / pi * norm(array_wt(:, j) * h), normalize(array_wt(:, j) * h)) * array_Rt(:, 3 * j - 2:3 * j);
    end

    for j = 1:num_links
        id1 = array_links(1, j);
        id2 = array_links(2, j);
        array_rt(1:3, j) = Rodrigues(180 / pi * norm(array_wt(:, id1) * h), normalize(array_wt(:, id1) * h)) * array_rt(1:3, j);
        array_rt(4:6, j) = Rodrigues(180 / pi * norm(array_wt(:, id2) * h), normalize(array_wt(:, id2) * h)) * array_rt(4:6, j);
    end

    for j = 1:num_bones
        x(3 * i - 2:3 * i, j) = array_xt(:, j);
        r(3 * i - 2:3 * i, 3 * j - 2:3 * j) = array_Rt(:, 3 * j - 2:3 * j);
    end
end