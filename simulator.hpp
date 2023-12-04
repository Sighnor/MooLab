#ifndef MOOLAB_SIMULATOR
#define MOOLAB_SIMULATOR

#include "motion.hpp"

typedef array2d<double> matrix;

static inline void print(matrix m1)
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

static inline matrix operator + (matrix m1, matrix m2)
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

static inline matrix operator - (matrix m1)
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

static inline matrix operator - (matrix m1, matrix m2)
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

static inline matrix operator * (double k, matrix m1)
{
    matrix m(m1.rows, m1.cols);
    for(int i = 0; i < m1.rows; i++)
    {
        for(int j = 0; j < m1.cols; j++)
        {
            m(i, j) = k * m1(i, j);
        }
    }
    return m;
}

static inline matrix operator * (matrix m1, matrix m2)
{
    assert(m1.cols == m2.rows);
    matrix m(m1.rows, m2.cols);
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

static inline matrix transpose_matrix(matrix m1)
{
    matrix m2(m1.rows, m1.cols);
    for(int i = 0; i < m2.rows; i++)
    {
        for(int j = 0; j < m2.cols; j++)
        {
            m2(i, j) = m1(j, i);
        }
    }
    return m2;
}

static inline void set_matrix(matrix &m1, int r1, int r2, int c1, int c2, matrix m2)
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

static inline void set_vec(vec3 &vec, matrix m1, int r1, int r2, int c1, int c2)
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
    for(int i = 0; i < m1.rows; i++)
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
    double radius;
    double length;
};


struct Simulator
{
    int bones_size;
    int links_size;

    array1d<int> bone_ids;
    array1d<int> bone_parents;
    array1d<vec3> bone_local_positions;
    array1d<quat> bone_local_rotations;

    array1d<double> bone_masses;
    array1d<mat3> bone_inertias;
    array1d<link> bone_links;
    array1d<shape> bone_shapes;

    array1d<vec3> bone_vels;
    array1d<vec3> bone_avels;
    array1d<vec3> bone_anim_positions;
    array1d<quat> bone_anim_rotations;

    array1d<vec3> bone_forces;
    array1d<vec3> bone_torgues;

    void bind_simulator(const BVH_Motion &motion);
    void simulate_gravity(double g);
    void simulate_damp(double damp);
    void simulate_contact();
    void simulate(double dt);
    void batch_forward_kinematics_full();
};

void Simulator::bind_simulator(const BVH_Motion &motion)
{
    bones_size = motion.nbones();
    bone_ids = motion.bone_ids;
    bone_parents = motion.bone_parents;
    bone_local_positions = motion.bone_local_positions(0);
    bone_local_rotations = motion.bone_local_rotations(0);

    bone_masses.resize(bones_size);
    bone_inertias.resize(bones_size);
    for(int i = 0; i < 1; i++)
    {
        bone_links.resize(bones_size);
    }
    bone_shapes.resize(bones_size);

    bone_vels.resize(bones_size);
    bone_avels.resize(bones_size);
    bone_anim_positions.resize(bones_size);
    bone_anim_rotations.resize(bones_size);

    bone_forces.resize(bones_size);
    bone_torgues.resize(bones_size);
}

void Simulator::simulate_gravity(double g)
{
    for(int i = 0; i < bones_size; i++)
    {
        bone_forces(i) = bone_forces(i) + vec3(0.f, 9.8 * (-bone_masses(i)), 0.f);
    }
}

void Simulator::simulate_damp(double damp)
{
    for(int i = 0; i < bones_size; i++)
    {
        bone_forces(i) = bone_forces(i) - damp * bone_masses(i) * bone_vels(i);
        bone_torgues(i) = bone_torgues(i) - damp * bone_masses(i) * bone_avels(i);
    }
}

void Simulator::simulate_contact()
{
    //碰撞检测
    for(int i = 0; i < bones_size; i++)
    {
        double eps = 1e-2;
        vec3 p1 = bone_anim_positions(i) + bone_anim_rotations(i) * vec3(0, 0.5 * bone_shapes(i).length - bone_shapes(i).radius, 0);
        vec3 p2 = bone_anim_positions(i) - bone_anim_rotations(i) * vec3(0, 0.5 * bone_shapes(i).length - bone_shapes(i).radius, 0);
        //刚体与刚体检测
        for(int j = i + 1; j < bones_size; j++)
        {
            vec3 q1 = bone_anim_positions(j) + bone_anim_rotations(j) * vec3(0, 0.5 * bone_shapes(j).length - bone_shapes(j).radius, 0);
            vec3 q2 = bone_anim_positions(j) - bone_anim_rotations(j) * vec3(0, 0.5 * bone_shapes(j).length - bone_shapes(j).radius, 0);

            double lambda1, lambda2;
            double distance = distance_between_line_segments(p1, p2, q1, q2, lambda1, lambda2);
            if(distance < bone_shapes(i).radius + bone_shapes(j).radius - eps)
            {
                float kp = 10000;
                float kd = 0.5;

                vec3 p = p1 + lambda1 * (p2 - p1);
                vec3 q = q1 + lambda2 * (q2 - q1);
                vec3 vp = (bone_vels(i) + cross(bone_avels(i), p - bone_anim_positions(i)));
                vec3 vq = (bone_vels(j) + cross(bone_avels(j), q - bone_anim_positions(j)));

                double length_Fp = kp * (bone_shapes(i).radius + bone_shapes(j).radius - eps - distance) + kd * dot(vp, q - p) / std::max(length(q - p), 1e-8f);
                double length_Fq = kp * (bone_shapes(i).radius + bone_shapes(j).radius - eps - distance) + kd * dot(vq, p - q) / std::max(length(p - q), 1e-8f);

                bone_forces(i) = bone_forces(i) + length_Fp * normalize(p - q);
                bone_torgues(i) = bone_torgues(i) + length_Fp * cross(p - bone_anim_positions(i), normalize(p - q));
                bone_forces(j) = bone_forces(j) + length_Fq * normalize(q - p);
                bone_torgues(j) = bone_torgues(j) + length_Fq * cross(q - bone_anim_positions(j), normalize(q - p));
            }
                
        }
        //地面碰撞检测
        double lambda;
        double distance = distance_between_line_segment_plane(p1, p2, 0, 1, 0, 1.1, lambda);
        if(distance < bone_shapes(i).radius - eps)
        {
            double kp = 10000;
            double kd = 100;

            vec3 p = p1 + lambda * (p2 - p1);
            vec3 v = (bone_vels(i) + cross(bone_avels(i), p - bone_anim_positions(i)));

            double length_F = kp * (bone_shapes(i).radius - eps - distance) + kd * dot(v, vec3(0, -1, 0));

            bone_forces(i) = bone_forces(i) + length_F * vec3(0, 1, 0);
            bone_torgues(i) = bone_torgues(i) + length_F * cross(p - bone_anim_positions(i), vec3(0, 1, 0));
        }  
    }
}

void Simulator::simulate(double h = 0.0166667)
{
    array1d<mat3> inertias(bones_size);

    matrix F(6 * bones_size, 1);
    matrix vt(6 * bones_size, 1);
    matrix M(6 * bones_size, 6 * bones_size);
    matrix J(3 * links_size, 6 * bones_size);
    matrix delta(3 * links_size, 1);

    for(int i = 0; i < bones_size; i++)
    {
        //计算实时转动惯量矩阵
        inertias(i) = quat_to_Rodrigues(bone_anim_rotations(i)) * bone_inertias(i) * trans_mat(quat_to_Rodrigues(bone_anim_rotations(i)));
        //更新重力和转动惯量阻力
        bone_torgues(i) = bone_torgues(i) - cross(bone_avels(i), bone_inertias(i) * bone_avels(i));
        set_matrix(F, 6 * i + 0, 6 * i + 2, 0, 0, bone_forces(i));
        set_matrix(F, 6 * i + 3, 6 * i + 5, 0, 0, bone_torgues(i));
        //更新速度
        set_matrix(vt, 6 * i + 0, 6 * i + 2, 0, 0, bone_vels(i));
        set_matrix(vt, 6 * i + 3, 6 * i + 5, 0, 0, bone_avels(i));
        //更新M矩阵
        set_matrix(M, 6 * i + 0, 6 * i + 2, 6 * i + 0, 6 * i + 2, 1 / h * bone_masses(i) * eye3());
        set_matrix(M, 6 * i + 3, 6 * i + 5, 6 * i + 3, 6 * i + 5, 1 / h * bone_inertias(i));
    }

    //计算实时误差和约束矩阵
    for(int i = 0; i < links_size; i++)
    {
        bone_links(i).delta = bone_anim_positions(bone_links(i).a) + bone_links(i).ra - bone_anim_positions(bone_links(i).b) - bone_links(i).rb;
        set_matrix(J, 3 * i + 0, 3 * i + 2, 6 * bone_links(i).a + 0, 6 * bone_links(i).a + 2, eye3());
        set_matrix(J, 3 * i + 0, 3 * i + 2, 6 * bone_links(i).a + 3, 6 * bone_links(i).a + 5, -vec_to_cross_matrix(bone_links(i).ra));
        set_matrix(J, 3 * i + 0, 3 * i + 2, 6 * bone_links(i).b + 0, 6 * bone_links(i).b + 2, -eye3());
        set_matrix(J, 3 * i + 0, 3 * i + 2, 6 * bone_links(i).b + 3, 6 * bone_links(i).b + 5, vec_to_cross_matrix(bone_links(i).rb));
        set_matrix(delta, 3 * i + 0, 3 * i + 2, 0, 0, -1 / h * bone_links(i).delta);
    }
    //解算方程
    matrix A(M.rows + J.rows, M.cols + J.rows);
    matrix B(M.rows + J.rows, 1);
    A.zero();
    B.zero();
    set_matrix(A, 0, M.rows - 1, 0, M.cols - 1, M);
    set_matrix(A, M.rows, M.rows + J.rows - 1, 0, M.cols - 1, J);
    set_matrix(A, 0, M.rows - 1, M.cols, M.cols + J.rows - 1, - transpose_matrix(J));

    set_matrix(B, 0, M.rows - 1, 0, 0, M * vt + F);
    set_matrix(B, M.rows, M.rows + J.rows - 1, 0, 0, delta);

    matrix vt_1 = Elimination_method_solution_of_linear_equations(A, B);
    //更新运动数据
    for(int i = 0; i < links_size; i++)
    {
        set_vec(bone_vels(i), vt_1, 6 * i + 0, 6 * i + 2, 0, 0);
        set_vec(bone_avels(i), vt_1, 6 * i + 3, 6 * i + 5, 0, 0);
        bone_anim_positions(i) = vel_to_vec(bone_vels(i), h) + bone_anim_positions(i);
        bone_anim_rotations(i) = avel_to_quat(bone_avels(i), h) * bone_anim_rotations(i);
    }
    //更新约束数据
    for(int i = 0; i < links_size; i++)
    {
        bone_links(i).ra = avel_to_quat(bone_avels(bone_links(i).a), h) * bone_links(i).ra;
        bone_links(i).rb = avel_to_quat(bone_avels(bone_links(i).b), h) * bone_links(i).rb;
    }
    //力和力矩归零
    for(int i = 0; i < bones_size; i++)
    {
        bone_forces(i) = vec3(0.f);
        bone_torgues(i) = vec3(0.f);
    }
}

void Simulator::batch_forward_kinematics_full()
{
    bone_anim_positions.zero();
    bone_anim_rotations.zero();
    
    for(int i = 0; i < bones_size; i++)
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

#endif