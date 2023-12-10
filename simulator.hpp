#ifndef MOOLAB_SIMULATOR
#define MOOLAB_SIMULATOR

#include "motion.hpp"

typedef array2d<double> matrix;

static inline void print(const matrix &m1)
{
    for(int i = 0; i < m1.rows; i++)
    {
        for(int j = 0; j < m1.cols; j++)
        {
            std::cout << m1(i, j) << " ";
        }
        std::cout << std::endl;
    }
}

static inline matrix operator + (const matrix &m1, const matrix &m2)
{
    assert(m1.rows == m2.rows && m1.cols == m2.cols);
    matrix m(m1.rows, m1.cols);
    for(int i = 0; i < m.rows; i++)
    {
        for(int j = 0; j < m.cols; j++)
        {
            m(i, j) = m1(i, j) + m2(i, j);
        }
    }
    return m;
}

static inline matrix operator - (const matrix &m1)
{
    matrix m(m1.rows, m1.cols);
    for(int i = 0; i < m.rows; i++)
    {
        for(int j = 0; j < m.cols; j++)
        {
            m(i, j) = - m1(i, j);
        }
    }
    return m;
}

static inline matrix operator - (const matrix &m1, const matrix &m2)
{
    assert(m1.rows == m2.rows && m1.cols == m2.cols);
    matrix m(m1.rows, m1.cols);
    for(int i = 0; i < m.rows; i++)
    {
        for(int j = 0; j < m.cols; j++)
        {
            m(i, j) = m1(i, j) - m2(i, j);
        }
    }
    return m;
}

static inline matrix operator * (const matrix &m1, const matrix &m2)
{
    assert(m1.cols == m2.rows);
    matrix m(m1.rows, m2.cols);
    // 必须归零，否则迭代中没有初值！
    m.zero();
    for(int i = 0; i < m1.rows; i++)
    {
        for(int j = 0; j < m2.cols; j++)
        {
            for(int k = 0; k < m1.cols; k++)
            {
                m(i, j) = m(i, j) + m1(i, k) * m2(k, j);
            }
        }
    }
    return m;
}

static inline matrix transpose_matrix(const matrix &m1)
{
    matrix m(m1.cols, m1.rows);
    for(int i = 0; i < m.rows; i++)
    {
        for(int j = 0; j < m.cols; j++)
        {
            m(i, j) = m1(j, i);
        }
    }
    return m;
}

static inline void set_matrix(matrix &m1, int r1, int r2, int c1, int c2, const matrix &m2)
{
    assert(r2 >= r1 && r1 >= 0 && r2 < m1.rows && c2 >= c1 && c1 >= 0 && c2 < m1.cols && r2 - r1 + 1 == m2.rows && c2 - c1 + 1 == m2.cols);
    for(int i = 0; i < m2.rows; i++)
    {
        for(int j = 0; j < m2.cols; j++)
        {
            m1(r1 + i, c1 + j) = m2(i, j);
        }
    }
}

static inline void set_matrix(matrix &m1, int r1, int r2, int c1, int c2, vec3 vec)
{
    assert(r2 >= r1 && r1 >= 0 && r2 < m1.rows && c2 >= c1 && c1 >= 0 && c2 < m1.cols && r2 - r1 + 1 == 3 && c2 - c1 + 1 == 1);
    m1(r1 + 0, c1) = vec.x;
    m1(r1 + 1, c1) = vec.y;
    m1(r1 + 2, c1) = vec.z;
}

static inline void set_matrix(matrix &m1, int r1, int r2, int c1, int c2, mat3 mat)
{
    assert(r2 >= r1 && r1 >= 0 && r2 < m1.rows && c2 >= c1 && c1 >= 0 && c2 < m1.cols && r2 - r1 + 1 == 3 && c2 - c1 + 1 == 3);
    m1(r1 + 0, c1 + 0) = mat.X.x;
    m1(r1 + 0, c1 + 1) = mat.Y.x;
    m1(r1 + 0, c1 + 2) = mat.Z.x;
    m1(r1 + 1, c1 + 0) = mat.X.y;
    m1(r1 + 1, c1 + 1) = mat.Y.y;
    m1(r1 + 1, c1 + 2) = mat.Z.y;
    m1(r1 + 2, c1 + 0) = mat.X.z;
    m1(r1 + 2, c1 + 1) = mat.Y.z;
    m1(r1 + 2, c1 + 2) = mat.Z.z;
}

static inline void set_vec(vec3 &vec, const matrix &m1, int r1, int r2, int c1, int c2)
{
    assert(r2 >= r1 && r1 >= 0 && r2 < m1.rows && c2 >= c1 && c1 >= 0 && c2 < m1.cols && r2 - r1 + 1 == 3 && c2 - c1 + 1 == 1);
    vec.x = m1(r1 + 0, c1);
    vec.y = m1(r1 + 1, c1);
    vec.z = m1(r1 + 2, c1);
}

static inline matrix concatenate_matrixs(const matrix &m1, const matrix &m2)
{
    assert(m1.rows == m2.rows);
    matrix m(m1.rows, m1.cols + m2.cols);
    for(int i = 0; i < m.rows; i++)
    {
        for(int j = 0; j < m1.cols; j++)
        {
            m(i, j) = m1(i, j);
        }
        for(int j = 0; j < m2.cols; j++)
        {
            m(i, m1.cols + j) = m2(i, j);
        }
    }
    return m;
}

static inline matrix Elimination_method_solution_of_linear_equations(const matrix &A, const matrix &y)
{
    assert(A.rows == A.cols && A.rows == y.rows && y.cols == 1);
    matrix B = concatenate_matrixs(A, y);
    matrix x(y.rows, y.cols);
    //对A除了最后一列的所有列执行消去(最后一列无需消去)
    for(int k = 0; k < B.cols - 2; k++)
    {
        //对角线下方每一个元素
        for(int i = k + 1; i < B.rows; i++)
        {
            //计算系数
            double m = B(i, k) / B(k, k);
            //行变换消去
            for(int j = 0; j < B.cols; j++)
            {
                B(i, j) = B(i, j) - m * B(k, j);
            }
        }
    }
        
    //从最后一行开始计算
    for(int i = B.rows - 1; i >= 0; i--)
    {
        //等式右方的值
        double temp = B(i, B.cols - 1);
        //回代减去其它x贡献
        for(int j = i + 1; j < B.cols - 1; j++)
        {
            temp = temp - B(i, j) * x(j, 0); 
        }
        //求值
        x(i, 0) = temp / B(i, i);
    }
    return x;
}

double distance_between_point_plane(vec3 p, double A, double B, double C, double D)
{
    return abs(A * p.x + B * p.y + C * p.z + D) / sqrt(A * A + B * B + C * C);
}

double distance_between_line_segment_plane(vec3 p1, vec3 p2, double A, double B, double C, double D, double &lambda)
{
    vec3 p(0.f, 0.f, 0.f);
    if(A != 0)
    {
        p = vec3(-D / A, 0, 0);
    }
    else if(B != 0)
    {
        p = vec3(0, -D / B, 0);
    }
    else if(C != 0)
    {
        p = vec3(0, 0, -D / C);
    }
    double dotp1pn = dot(p1 - p, vec3(A, B, C));
    double dotp2pn = dot(p2 - p, vec3(A, B, C));
    double distancep1plane = sgn(dotp1pn) * distance_between_point_plane(p1, A, B, C, D);
    double distancep2plane = sgn(dotp2pn) * distance_between_point_plane(p2, A, B, C, D);
   
    if(distancep1plane <= distancep2plane)
    {
        lambda = 0;
        return distancep1plane;
    }
    else
    {
        lambda = 1;
        return distancep2plane;
    }
    return 0.f;
}

double distance_between_point_line_segment(vec3 p0, vec3 p1, vec3 p2, double &lambda)
{
    vec3 s = p2 - p1;
    double den = dot(s, s);
    lambda = dot((p0 - p1), s) / std::max(den, 1e-8);
    if(lambda >= 0 && lambda <= 1 && den >= 1e-8)
    {
        return length(p0 - (p1 + lambda * s));
    }
    else
    {
        double distance_p0p1 = length(p0 - p1);
        double distance_p0p2 = length(p0 - p2);
        if(distance_p0p1 <= distance_p0p2)
        {
            lambda = 0;
            return distance_p0p1;
        }
        else
        {
            lambda = 1;
            return distance_p0p2;
        }    
    }
    return 0.f;
}

double distance_between_line_segments(vec3 p1, vec3 p2, vec3 q1, vec3 q2, double &lambda1, double &lambda2)
{
    vec3 s1 = p2 - p1;
    vec3 s2 = q2 - q1;
    double den = (dot(s1, s1) * dot(s2, s2) - dot(s1, s2) * dot(s1, s2));
    lambda1 = (dot(s1, s2) * dot((p1 - q1), s2) - dot(s2, s2) * dot(p1 - q1, s1)) / std::max(den, 1e-8);
    lambda2 = -(dot(s1, s2) * dot((p1 - q1), s1) - dot(s1, s1) * dot(p1 - q1, s2)) / std::max(den, 1e-8);
    if(lambda1 >= 0 && lambda1 <= 1 && lambda2 >= 0 && lambda2 <= 1 && den >= 1e-8)
    {
        return length(p1 + lambda1 * s1 - q1 - lambda2 * s2);
    }
    else
    {
        double lambda_p1l2, lambda_p2l2, lambda_q1l1, lambda_q2l1;
        double distance_p1l2 = distance_between_point_line_segment(p1, q1, q2, lambda_p1l2);
        double distance_p2l2 = distance_between_point_line_segment(p2, q1, q2, lambda_p2l2);
        double distance_q1l1 = distance_between_point_line_segment(q1, p1, p2, lambda_q1l1);
        double distance_q2l1 = distance_between_point_line_segment(q2, p1, p2, lambda_q2l1);
        if(distance_p1l2 <= distance_p2l2 && distance_p1l2 <= distance_q1l1 && distance_p1l2 <= distance_q2l1)
        {
            lambda1 = 0;
            lambda2 = lambda_p1l2;
            return distance_p1l2;
        }
        else if(distance_p2l2 <= distance_q1l1 && distance_p2l2 <= distance_q2l1)
        {
            lambda1 = 1;
            lambda2 = lambda_p2l2;
            return distance_p2l2;
        }
        else if(distance_q1l1 <= distance_q2l1)
        {
            lambda1 = lambda_q1l1;
            lambda2 = 0;
            return distance_q1l1;
        }
        else
        {
            lambda1 = lambda_q2l1;
            lambda2 = 1;
            return distance_q2l1;
        }  
    }
    return 0.f;
}

struct link
{
    int a;
    int b;
    vec3 ra;
    vec3 rb;
    vec3 delta;
};

struct shape
{
    vec3 pos;
    quat rot;
    double radius;
    double length;
};


struct Simulator
{
    array1d<int> bone_ids;
    array1d<int> bone_parents;
    array1d<vec3> bone_local_positions;
    array1d<quat> bone_local_rotations;
    array1d<vec3> bone_anim_positions;
    array1d<quat> bone_anim_rotations;
    array1d<quat> bone_quat_compensates;

    array1d<double> bone_masses;
    array1d<mat3> bone_inertias;
    array1d<link> bone_links;
    array1d<shape> bone_shapes;

    array1d<vec3> bone_vels;
    array1d<vec3> bone_avels;

    array1d<vec3> bone_forces;
    array1d<vec3> bone_torgues;

    void simulate_gravity(double g);
    void simulate_damp(double damp);
    void simulate_contact(double ground_height, double eps);
    void simulate(double dt);
    void batch_forward_kinematics_full();
};

void Simulator::simulate_gravity(double g)
{
    for(int i = 0; i < bone_masses.size; i++)
    {
        bone_forces(i) = bone_forces(i) + vec3(0.f, -9.8f * bone_masses(i), 0.f);
    }
}

void Simulator::simulate_damp(double damp)
{
    for(int i = 0; i < bone_masses.size; i++)
    {
        bone_forces(i) = bone_forces(i) - damp * bone_masses(i) * bone_vels(i);
        bone_torgues(i) = bone_torgues(i) - damp * bone_masses(i) * bone_avels(i);
    }
}

void Simulator::simulate_contact(double ground_height = 0, double eps = 1e-2)
{
    //碰撞检测
    for(int i = 0; i < bone_masses.size; i++)
    {
        vec3 p1 = bone_shapes(i).pos + bone_shapes(i).rot * vec3(0.f, 0.5f * bone_shapes(i).length - bone_shapes(i).radius, 0.f);
        vec3 p2 = bone_shapes(i).pos - bone_shapes(i).rot * vec3(0.f, 0.5f * bone_shapes(i).length - bone_shapes(i).radius, 0.f);
        //刚体与刚体检测
        for(int j = i + 1; j < bone_masses.size; j++)
        {
            vec3 q1 = bone_shapes(j).pos + bone_shapes(j).rot * vec3(0.f, 0.5f * bone_shapes(j).length - bone_shapes(j).radius, 0.f);
            vec3 q2 = bone_shapes(j).pos - bone_shapes(j).rot * vec3(0.f, 0.5f * bone_shapes(j).length - bone_shapes(j).radius, 0.f);

            double lambda1, lambda2;
            double distance = distance_between_line_segments(p1, p2, q1, q2, lambda1, lambda2);
            if(distance < bone_shapes(i).radius + bone_shapes(j).radius - eps)
            {
                float kp = 10000;
                float kd = 0.5;

                vec3 p = p1 + lambda1 * (p2 - p1);
                vec3 q = q1 + lambda2 * (q2 - q1);
                vec3 vp = (bone_vels(i) + cross(bone_avels(i), p - bone_shapes(i).pos));
                vec3 vq = (bone_vels(j) + cross(bone_avels(j), q - bone_shapes(j).pos));

                double length_Fp = kp * (bone_shapes(i).radius + bone_shapes(j).radius - eps - distance) + kd * dot(vp, q - p) / std::max(length(q - p), 1e-8f);
                double length_Fq = kp * (bone_shapes(i).radius + bone_shapes(j).radius - eps - distance) + kd * dot(vq, p - q) / std::max(length(p - q), 1e-8f);

                bone_forces(i) = bone_forces(i) + length_Fp * normalize(p - q);
                bone_torgues(i) = bone_torgues(i) + length_Fp * cross(p - bone_shapes(i).pos, normalize(p - q));
                bone_forces(j) = bone_forces(j) + length_Fq * normalize(q - p);
                bone_torgues(j) = bone_torgues(j) + length_Fq * cross(q - bone_shapes(j).pos, normalize(q - p));
            }
                
        }
        //地面碰撞检测
        double lambda;
        double distance = distance_between_line_segment_plane(p1, p2, 0, 1, 0, -ground_height, lambda);
        if(distance < bone_shapes(i).radius - eps)
        {
            double kp = 10000;
            double kd = 100;

            vec3 p = p1 + lambda * (p2 - p1);
            vec3 v = (bone_vels(i) + cross(bone_avels(i), p - bone_shapes(i).pos));

            double length_F = kp * (bone_shapes(i).radius - eps - distance) + kd * dot(v, vec3(0, -1, 0));

            bone_forces(i) = bone_forces(i) + length_F * vec3(0, 1, 0);
            bone_torgues(i) = bone_torgues(i) + length_F * cross(p - bone_shapes(i).pos, vec3(0, 1, 0));
        }  
    }
}

void Simulator::simulate(double h = 0.0166667)
{
    array1d<mat3> inertias(bone_masses.size);

    matrix F(6 * bone_masses.size, 1);
    matrix vt(6 * bone_masses.size, 1);
    matrix M(6 * bone_masses.size, 6 * bone_masses.size);
    matrix J(3 * bone_links.size, 6 * bone_masses.size);
    matrix delta(3 * bone_links.size, 1);
    F.zero();
    vt.zero();
    M.zero();
    J.zero();
    delta.zero();

    for(int i = 0; i < bone_masses.size; i++)
    {
        //计算实时转动惯量矩阵
        inertias(i) = quat_to_Rodrigues(bone_shapes(i).rot) * bone_inertias(i) * trans_mat(quat_to_Rodrigues(bone_shapes(i).rot));
        //更新转动惯量阻力
        bone_torgues(i) = bone_torgues(i) - cross(bone_avels(i), inertias(i) * bone_avels(i));
        set_matrix(F, 6 * i + 0, 6 * i + 2, 0, 0, bone_forces(i));
        set_matrix(F, 6 * i + 3, 6 * i + 5, 0, 0, bone_torgues(i));
        //更新速度
        set_matrix(vt, 6 * i + 0, 6 * i + 2, 0, 0, bone_vels(i));
        set_matrix(vt, 6 * i + 3, 6 * i + 5, 0, 0, bone_avels(i));
        //更新M矩阵
        set_matrix(M, 6 * i + 0, 6 * i + 2, 6 * i + 0, 6 * i + 2, 1 / h * bone_masses(i) * eye3());
        set_matrix(M, 6 * i + 3, 6 * i + 5, 6 * i + 3, 6 * i + 5, 1 / h * inertias(i));
    }

    //计算实时误差和约束矩阵
    for(int i = 0; i < bone_links.size; i++)
    {
        int id1 = bone_links(i).a;
        int id2 = bone_links(i).b;
        bone_links(i).delta = bone_shapes(id1).pos + bone_links(i).ra - bone_shapes(id2).pos - bone_links(i).rb;
        // print(bone_links(i).delta);
        set_matrix(J, 3 * i + 0, 3 * i + 2, 6 * id1 + 0, 6 * id1 + 2, eye3());
        set_matrix(J, 3 * i + 0, 3 * i + 2, 6 * id1 + 3, 6 * id1 + 5, -vec_to_cross_matrix(bone_links(i).ra));
        set_matrix(J, 3 * i + 0, 3 * i + 2, 6 * id2 + 0, 6 * id2 + 2, -eye3());
        set_matrix(J, 3 * i + 0, 3 * i + 2, 6 * id2 + 3, 6 * id2 + 5, vec_to_cross_matrix(bone_links(i).rb));
        set_matrix(delta, 3 * i + 0, 3 * i + 2, 0, 0, -1 / h * bone_links(i).delta);
    }
    //解算方程
    matrix A(M.rows + J.rows, M.cols + J.rows);
    matrix B(M.rows + J.rows, 1);
    matrix M_vt_F = M * vt + F;
    matrix neg_trans_J = -transpose_matrix(J);
    A.zero();
    B.zero();

    set_matrix(A, 0, M.rows - 1, 0, M.cols - 1, M);
    set_matrix(A, M.rows, M.rows + J.rows - 1, 0, M.cols - 1, J);
    set_matrix(A, 0, M.rows - 1, M.cols, M.cols + J.rows - 1, neg_trans_J);

    set_matrix(B, 0, M.rows - 1, 0, 0, M_vt_F);
    set_matrix(B, M.rows, M.rows + J.rows - 1, 0, 0, delta);

    matrix vt_1 = Elimination_method_solution_of_linear_equations(A, B);
    //更新运动数据
    for(int i = 0; i < bone_masses.size; i++)
    {
        set_vec(bone_vels(i), vt_1, 6 * i + 0, 6 * i + 2, 0, 0);
        set_vec(bone_avels(i), vt_1, 6 * i + 3, 6 * i + 5, 0, 0);
        // print(bone_vels(i));
        bone_shapes(i).pos = vel_to_vec(bone_vels(i), h) + bone_shapes(i).pos;
        bone_shapes(i).rot = avel_to_quat(bone_avels(i), h) * bone_shapes(i).rot;
    }
    //更新约束数据
    for(int i = 0; i < bone_links.size; i++)
    {
        int id1 = bone_links(i).a;
        int id2 = bone_links(i).b;
        bone_links(i).ra = avel_to_quat(bone_avels(id1), h) * bone_links(i).ra;
        bone_links(i).rb = avel_to_quat(bone_avels(id2), h) * bone_links(i).rb;
    }
    //力和力矩归零
    for(int i = 0; i < bone_masses.size; i++)
    {
        bone_forces(i) = vec3(0.f);
        bone_torgues(i) = vec3(0.f);
    }
}

void Simulator::batch_forward_kinematics_full()
{
    bone_anim_positions.zero();
    bone_anim_rotations.zero();
    
    for(int i = 0; i < bone_ids.size; i++)
    {
        // Assumes bones are always sorted from root onwards
        int parent_id = bone_parents(i);

        if(parent_id == -1)
        {
            bone_anim_rotations(i) = bone_local_rotations(i);
            bone_anim_positions(i) = bone_local_positions(i);
        }
        else
        {
            bone_anim_rotations(i) = 
                bone_anim_rotations(parent_id) * bone_local_rotations(i);
            bone_anim_positions(i) = 
                bone_anim_positions(parent_id) + 
                bone_anim_rotations(parent_id) * bone_local_positions(i);
        }
    }
}

void bind_simulator(Simulator &simulator, const BVH_Motion &motion, int frame_id, double radius)
{
    simulator.bone_ids = motion.bone_ids;
    simulator.bone_parents = motion.bone_parents;
    simulator.bone_local_positions = motion.bone_local_positions(frame_id);
    simulator.bone_local_rotations = motion.bone_local_rotations(frame_id);
    simulator.bone_anim_positions.resize(motion.nbones());
    simulator.bone_anim_rotations.resize(motion.nbones());
    simulator.bone_quat_compensates.resize(motion.nbones());

    simulator.bone_masses.resize(motion.nbones() - 1);
    simulator.bone_inertias.resize(motion.nbones() - 1);
    simulator.bone_links.resize(motion.nbones() - 2);
    simulator.bone_shapes.resize(motion.nbones() - 1);

    simulator.bone_vels.resize(motion.nbones() - 1);
    simulator.bone_vels.zero();
    simulator.bone_avels.resize(motion.nbones() - 1);
    simulator.bone_avels.zero();

    simulator.bone_forces.resize(motion.nbones() - 1);
    simulator.bone_forces.zero();
    simulator.bone_torgues.resize(motion.nbones() - 1);
    simulator.bone_torgues.zero();

    for(int i = 0; i < simulator.bone_masses.size; i++)
    {
        simulator.bone_masses(i) = 1;
        simulator.bone_inertias(i) = eye3();
    }

    simulator.batch_forward_kinematics_full();

    for(int i = 2; i < simulator.bone_ids.size; i++)
    {
        if(i != 2 && i != 6 && i != 10 && i != 13 && i != 15 && i != 19)
        {
            int parent_id = simulator.bone_parents(i);
            simulator.bone_quat_compensates(parent_id - 1) = vec_to_quat(vec3(0.f, 1.f, 0.f), simulator.bone_local_positions(i));
            simulator.bone_shapes(parent_id - 1).pos = 0.5f * (simulator.bone_anim_positions(i) + (simulator.bone_anim_positions(parent_id)));
            simulator.bone_shapes(parent_id - 1).rot = simulator.bone_anim_rotations(parent_id) * simulator.bone_quat_compensates(parent_id - 1);
            simulator.bone_shapes(parent_id - 1).radius = radius;
            simulator.bone_shapes(parent_id - 1).length = length(simulator.bone_anim_positions(i) - (simulator.bone_anim_positions(parent_id)));
        }
    }

    simulator.bone_quat_compensates(1 - 1) = vec_to_quat(vec3(0.f, 1.f, 0.f), simulator.bone_local_positions(1));
    simulator.bone_shapes(1 - 1).pos = simulator.bone_anim_positions(1);
    simulator.bone_shapes(1 - 1).rot = simulator.bone_anim_rotations(1) * simulator.bone_quat_compensates(1 - 1);
    simulator.bone_shapes(1 - 1).radius = minf_3(
                                            length(simulator.bone_shapes(1 - 1).pos - simulator.bone_anim_positions(2)), 
                                            length(simulator.bone_shapes(1 - 1).pos - simulator.bone_anim_positions(6)), 
                                            length(simulator.bone_shapes(1 - 1).pos - simulator.bone_anim_positions(10)));
    simulator.bone_shapes(1 - 1).length = 2 * simulator.bone_shapes(1 - 1).radius;

    simulator.bone_quat_compensates(12 - 1) = vec_to_quat(vec3(0.f, 1.f, 0.f), simulator.bone_local_positions(12));
    simulator.bone_shapes(12 - 1).pos = 0.5 * (simulator.bone_anim_positions(12) + simulator.bone_anim_positions(13));
    simulator.bone_shapes(12 - 1).rot = simulator.bone_anim_rotations(12) * simulator.bone_quat_compensates(12 - 1);
    simulator.bone_shapes(12 - 1).radius = std::min(
                                                std::min(
                                                    length(simulator.bone_shapes(12 - 1).pos - simulator.bone_anim_positions(12)), 
                                                    length(simulator.bone_shapes(12 - 1).pos - simulator.bone_anim_positions(13))), 
                                                std::min(
                                                    length(simulator.bone_shapes(12 - 1).pos - simulator.bone_anim_positions(15)), 
                                                    length(simulator.bone_shapes(12 - 1).pos - simulator.bone_anim_positions(19))));
    simulator.bone_shapes(12 - 1).length = 2 * simulator.bone_shapes(12 - 1).radius;

    simulator.bone_quat_compensates(5 - 1) = vec_to_quat(vec3(0.f, 1.f, 0.f), simulator.bone_local_positions(5));
    simulator.bone_shapes(5 - 1).pos = simulator.bone_anim_positions(5) + simulator.bone_anim_rotations(5) * vec3(radius, 0.f, 0.f);
    simulator.bone_shapes(5 - 1).rot = simulator.bone_anim_rotations(5) * simulator.bone_quat_compensates(5 - 1);
    simulator.bone_shapes(5 - 1).radius = radius;
    simulator.bone_shapes(5 - 1).length = 2 * radius;

    simulator.bone_quat_compensates(9 - 1) = vec_to_quat(vec3(0.f, 1.f, 0.f), simulator.bone_local_positions(9));
    simulator.bone_shapes(9 - 1).pos = simulator.bone_anim_positions(9) + simulator.bone_anim_rotations(9) * vec3(radius, 0.f, 0.f);
    simulator.bone_shapes(9 - 1).rot = simulator.bone_anim_rotations(9) * simulator.bone_quat_compensates(9 - 1);
    simulator.bone_shapes(9 - 1).radius = radius;
    simulator.bone_shapes(9 - 1).length = 2 * radius;

    simulator.bone_quat_compensates(14 - 1) = vec_to_quat(vec3(0.f, 1.f, 0.f), simulator.bone_local_positions(14));
    simulator.bone_shapes(14 - 1).pos = simulator.bone_anim_positions(14) + simulator.bone_anim_rotations(14) * vec3(radius, 0.f, 0.f);
    simulator.bone_shapes(14 - 1).rot = simulator.bone_anim_rotations(14) * simulator.bone_quat_compensates(14 - 1);
    simulator.bone_shapes(14 - 1).radius = radius;
    simulator.bone_shapes(14 - 1).length = 2 * radius;

    simulator.bone_quat_compensates(18 - 1) = vec_to_quat(vec3(0.f, 1.f, 0.f), simulator.bone_local_positions(18));
    simulator.bone_shapes(18 - 1).pos = simulator.bone_anim_positions(18) + simulator.bone_anim_rotations(18) * vec3(radius, 0.f, 0.f);
    simulator.bone_shapes(18 - 1).rot = simulator.bone_anim_rotations(18) * simulator.bone_quat_compensates(18 - 1);
    simulator.bone_shapes(18 - 1).radius = radius;
    simulator.bone_shapes(18 - 1).length = 2 * radius;

    simulator.bone_quat_compensates(22 - 1) = vec_to_quat(vec3(0.f, 1.f, 0.f), simulator.bone_local_positions(22));
    simulator.bone_shapes(22 - 1).pos = simulator.bone_anim_positions(22) + simulator.bone_anim_rotations(22) * vec3(radius, 0.f, 0.f);
    simulator.bone_shapes(22 - 1).rot = simulator.bone_anim_rotations(22) * simulator.bone_quat_compensates(22 - 1);
    simulator.bone_shapes(22 - 1).radius = radius;
    simulator.bone_shapes(22 - 1).length = 2 * radius;

    for(int i = 2; i < simulator.bone_ids.size; i++)
    {
        int parent_id = simulator.bone_parents(i);
        simulator.bone_links(i - 2).a = parent_id - 1;
        simulator.bone_links(i - 2).b = i - 1;
        simulator.bone_links(i - 2).ra = simulator.bone_anim_positions(i) - simulator.bone_shapes(parent_id - 1).pos;
        simulator.bone_links(i - 2).rb = simulator.bone_anim_positions(i) - simulator.bone_shapes(i - 1).pos;
        simulator.bone_links(i - 2).delta = 0;
    }

    // for(int i = 0; i < simulator.bone_masses.size; i++)
    // {
    //     print(simulator.bone_shapes(i).pos);
    // }
    // for(int i = 0; i < simulator.bone_masses.size; i++)
    // {
    //     print(simulator.bone_shapes(i).rot);
    // }
    // for(int i = 0; i < simulator.bone_masses.size; i++)
    // {
    //     std::cout << (simulator.bone_shapes(i).radius) << ", " << (simulator.bone_shapes(i).length) << std::endl;
    // }
    // for(int i = 0; i < simulator.bone_links.size; i++)
    // {
    //     std::cout << (simulator.bone_links(i).a) << ", " 
    //               << (simulator.bone_links(i).b) << std::endl;
    // }
    // for(int i = 0; i < simulator.bone_links.size; i++)
    // {
    //     std::cout << (simulator.bone_links(i).ra.x) << ", " 
    //               << (simulator.bone_links(i).ra.y) << ", " 
    //               << (simulator.bone_links(i).ra.z) << ", " 
    //               << (simulator.bone_links(i).rb.x) << ", " 
    //               << (simulator.bone_links(i).rb.y) << ", " 
    //               << (simulator.bone_links(i).rb.z) << std::endl;
    // }
}

void deform_character_anim_bones(const Simulator &simulator, slice1d<vec3> _bone_anim_positions, slice1d<quat> _bone_anim_rotations)
{
    _bone_anim_positions.zero();
    _bone_anim_rotations.zero();
    
    for(int i = 1; i < simulator.bone_ids.size; i++)
    {
        _bone_anim_rotations(i) = simulator.bone_shapes(i - 1).rot * inv_quat(simulator.bone_quat_compensates(i - 1));
    }

    _bone_anim_positions(1) = simulator.bone_shapes(0).pos;
    for(int i = 0; i < simulator.bone_links.size; i++)
    {
        int id1 = simulator.bone_links(i).a;
        int id2 = simulator.bone_links(i).b;
        _bone_anim_positions(id2 + 1) = (simulator.bone_shapes(id1).pos + simulator.bone_links(i).ra + simulator.bone_shapes(id2).pos + simulator.bone_links(i).rb) / 2.f;
    }
}

#endif