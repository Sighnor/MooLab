#ifndef MOOLAB_CONTROLLER
#define MOOLAB_CONTROLLER

#include "quat.hpp"
#include "vec.hpp"

struct Camera_Controller
{
    vec3 *pos;
    vec3 vel;
    vec3 acc;

    vec3 *dir;
    vec2 ang;
    vec2 ang_vel;
    vec2 ang_acc;

    void pos_pid_control(
            float p_const, 
            float i_const, 
            float d_const, 
            float k, 
            float dt, 
            vec3 input);

    void dir_pid_control(
            float p_const, 
            float i_const, 
            float d_const, 
            float k, 
            float dt, 
            vec3 input);

    void pos_gamepad_control(vec3 input, mat3 ori, float distance);
    void dir_gamepad_control(vec3 input, float k);
};

void Camera_Controller::pos_pid_control(
                    float p_const, 
                    float i_const, 
                    float d_const, 
                    float k, 
                    float dt, 
                    vec3 input)
{
    static vec3 pos_ek_0 = vec3(0.f);
    static vec3 pos_ek_1 = vec3(0.f);
    static vec3 pos_ek_2 = vec3(0.f);

    pos_ek_2 = pos_ek_1;
    pos_ek_1 = pos_ek_0;
    pos_ek_0 = input - *pos;

    acc = k * ((p_const * (pos_ek_0 - pos_ek_1) + i_const * pos_ek_0 - d_const * (pos_ek_0 - 2 * pos_ek_1 + pos_ek_2)));
    vel = vel + acc * dt;
    *pos = *pos + vel * dt;
}

void Camera_Controller::dir_pid_control(
                    float p_const, 
                    float i_const, 
                    float d_const, 
                    float k, 
                    float dt, 
                    vec3 input)
{
    static vec2 ang_ek_0 = vec2(0.f);
    static vec2 ang_ek_1 = vec2(0.f);
    static vec2 ang_ek_2 = vec2(0.f);

    vec2 target_ang = dir_to_sph(input);

    ang_ek_2 = ang_ek_1;
    ang_ek_1 = ang_ek_0;
    ang_ek_0 = target_ang - ang;

    ang_ek_0.x = circulate_float(ang_ek_0.x, - PI, PI);

    ang_acc = k * (p_const * (ang_ek_0 - ang_ek_1) + i_const * ang_ek_0 - d_const * (ang_ek_0 - 2 * ang_ek_1 + ang_ek_2));
    ang_vel = ang_vel + ang_acc * dt;
    ang = ang + ang_vel * dt;

    ang.x = circulate_float(ang.x, - PI, PI);
    ang.y = circulate_float(ang.y, 0.f, PI);

    *dir = sph_to_dir(ang);
}

void Camera_Controller::pos_gamepad_control(vec3 input = vec3(0.f), mat3 ori = eye3(), float distance = 3.f)
{
    // 默认方向为(0, 0, 1)，故默认向量为(0, 0 , -1)；
    // *pos = input + quat(rad_to_deg(ang.x), vec3(0.f, 1.f, 0.f)) * vec3(0.f, 0.f, -distance);
    *pos = input + ori * vec3(0.f, 0.f, distance);
}

void Camera_Controller::dir_gamepad_control(vec3 input, float k)
{
    ang.x = circulate_float(ang.x - k * input.x, - PI, PI);
    ang.y = clampf(ang.y + k * input.z, PI / 2, PI - 0.1f);
    *dir = sph_to_dir(ang);
}

void bind_controller(Camera_Controller &c, vec3 *pos, vec3 *dir)
{
    c.pos = pos;
    c.dir = dir;
    c.ang = dir_to_sph(*dir);
}

enum Bait
{
    stand = 0,
    walk  = 1,
    run   = 2,
    jump  = 3,
    squat = 4,
    roll  = 5
};

struct Character_Controller
{
    vec3 vel0;
    vec3 vel1;
    vec3 vel2;
    vec3 vel3;
    vec3 vel4;
    vec3 vel5;

    vec3 avel0;
    vec3 avel1;
    vec3 avel2;
    vec3 avel3;
    vec3 avel4;
    vec3 avel5;

    quat rotation0;

    Bait bait;

    void update(vec3 vel, vec3 avel, quat rotation);
};

void Character_Controller::update(vec3 vel, vec3 avel, quat rotation)
{
    vel0 = inv_quat(rotation) * rotation0 * vel5;
    avel0 = inv_quat(rotation) * rotation0 * avel5;
    // vel0 = vel5;
    // avel0 = avel5;
    vel5 = vel;
    avel5 = avel;
    rotation0 = rotation;

    vel1 = 0.8f * vel0 + 0.2f * vel5;
    vel2 = 0.6f * vel0 + 0.4f * vel5;
    vel3 = 0.4f * vel0 + 0.6f * vel5;
    vel4 = 0.2f * vel0 + 0.8f * vel5;

    avel1 = 0.8f * avel0 + 0.2f * avel5;
    avel2 = 0.6f * avel0 + 0.4f * avel5;
    avel3 = 0.4f * avel0 + 0.6f * avel5;
    avel4 = 0.2f * avel0 + 0.8f * avel5;
}

struct formation
{
    // 周期，记录每个队形点对应的周期内时间
    array1d<float> phases;
    // 记录每个队形点对应的位置
    array1d<vec3> poses;
    array1d<vec3> vels;
    array1d<vec3> accs;
    // 记录每个队形点对应的期望位置
    array1d<vec3> target_poses;
    array1d<vec3> poses_ek_0;
    array1d<vec3> poses_ek_1;
    array1d<vec3> poses_ek_2;
    // 对应的姿态，不过此处并未用到
    array1d<quat> rots;
    array1d<vec3> avels;
    array1d<vec3> aaccs;
    array1d<quat> target_rots;
    array1d<quat> rots_ek_0;
    array1d<quat> rots_ek_1;
    array1d<quat> rots_ek_2;

    void pid_control( 
            float p_const, 
            float i_const, 
            float d_const, 
            float k, 
            float dt);
    void formation_circle(
            vec3 c, 
            float r, 
            float cycle, 
            quat R, 
            float dt);
    void formation_triangle(
            vec3 c, 
            float length, 
            float cycle, 
            quat R, 
            float dt);
};

formation formation_init(int size, vec3 c, float r)
{
    formation f;

    f.phases.resize(size);

    f.poses.resize(size);
    for(int i = 0; i < size; i++)
    {
        f.phases(i) = i / float(size);
        f.poses(i) = c + vec3(r * sin(f.phases(i) * 2 * PI), 0.f, r * cos(f.phases(i) * 2 * PI));
    }
    f.vels.resize(size);
    f.vels.zero();
    f.accs.resize(size);
    f.accs.zero();
    f.target_poses.resize(size);
    f.target_poses.zero();
    f.poses_ek_0.resize(size);
    f.poses_ek_0.zero();
    f.poses_ek_1.resize(size);
    f.poses_ek_1.zero();
    f.poses_ek_2.resize(size);
    f.poses_ek_2.zero();

    f.rots.resize(size);
    for(int i = 0; i < size; i++)
    {
        f.phases(i) = i / float(size);
        f.rots(i) = quat(rad_to_deg(circulate_float(f.phases(i) * 2 * PI + PI / 2.f, - PI, PI)), vec3(0.f, 1.f, 0.f));
    }
    f.avels.resize(size);
    f.avels.zero();
    f.aaccs.resize(size);
    f.aaccs.zero();
    f.target_rots.resize(size);
    f.target_rots.zero();
    f.rots_ek_0.resize(size);
    f.rots_ek_1.resize(size);
    f.rots_ek_2.resize(size);
    for(int i = 0; i < size; i++)
    {
        f.rots_ek_0(i) = quat(1.f, 0.f, 0.f, 0.f);
        f.rots_ek_1(i) = quat(1.f, 0.f, 0.f, 0.f);
        f.rots_ek_2(i) = quat(1.f, 0.f, 0.f, 0.f);
    }

    return f;
}

void formation::pid_control( 
            float p_const, 
            float i_const, 
            float d_const, 
            float k, 
            float dt)
{
    for(int i = 0; i < poses.size; i++)
    {
        // 位置PID计算
        poses_ek_2(i) = poses_ek_1(i);
        poses_ek_1(i) = poses_ek_0(i);
        poses_ek_0(i) = target_poses(i) - poses(i);
        accs(i) = k * ((p_const * (poses_ek_0(i) - poses_ek_1(i)) + i_const * poses_ek_0(i) - d_const * (poses_ek_0(i) - 2 * poses_ek_1(i) + poses_ek_2(i))));
        vels(i) = vels(i) + accs(i) * dt;
        poses(i) = poses(i) + vels(i) * dt;
        // 姿态PID计算，此处并没用到
        rots_ek_2(i) = rots_ek_1(i);
        rots_ek_1(i) = rots_ek_0(i);
        rots_ek_0(i) = target_rots(i) * inv_quat(rots(i));
        aaccs(i) = 0.1 * k * ((p_const * quat_to_avel(rots_ek_1(i), rots_ek_0(i), 1) + 0 * i_const * quat_to_avel(quat(1.f, 0.f, 0.f, 0.f), rots_ek_0(i), 1) - 0 * d_const * quat_to_avel(rots_ek_1(i) * inv_quat(rots_ek_2(i)), (rots_ek_0(i) * inv_quat(rots_ek_1(i))), 1)));
        avels(i) = avels(i) + aaccs(i) * dt;
        rots(i) = avel_to_quat(avels(i), dt) * rots(i);
        rots(i) = target_rots(i);
    }
}

void formation::formation_circle( 
                    vec3 c, 
                    float r, 
                    float cycle, 
                    quat R, 
                    float dt)
{
    for(int i = 0; i < poses.size; i++)
    {
        phases(i) = circulate_float(phases(i) + dt / cycle, 0, 1);
        target_poses(i) = c + R * vec3(r * sin(phases(i) * 2 * PI), 0.f, r * cos(phases(i) * 2 * PI));
        target_rots(i) = R * quat(rad_to_deg(circulate_float(phases(i) * 2 * PI + PI / 2.f, - PI, PI)), vec3(0.f, 1.f, 0.f));
    }
}

void formation::formation_triangle(
                    vec3 c, 
                    float length, 
                    float cycle, 
                    quat R, 
                    float dt)
{
    for(int i = 0; i < poses.size; i++)
    {
        phases(i) = circulate_float(phases(i) + dt / cycle, 0, 1);
        if(phases(i) <= 2 / 9.f)
        {
            target_poses(i) = c + R * (length / sqrt(3) * (vec3(0, 0, 1) + phases(i) / (2 / 9.f) * vec3(sqrt(3), 0, - 1)));
            target_rots(i) = R * quat(120, vec3(0.f, 1.f, 0.f));
        }
        else if(phases(i) > 2 / 9.f && phases(i) <= 5 / 9.f)
        {
            target_poses(i) = c + R * (length / sqrt(3) * (vec3(sqrt(3), 0, 0) + (phases(i) - 2 / 9.f) / (1 / 3.f) * vec3(- 3 * sqrt(3) / 2.f, 0, - 3 / 2.f)));
            target_rots(i) = R * quat(-120, vec3(0.f, 1.f, 0.f));
        } 
        else if(phases(i) > 5 / 9.f && phases(i) <= 8 / 9.f)
        {
            target_poses(i) = c + R * (length / sqrt(3) * (vec3(- sqrt(3) / 2.f, 0, -3 / 2.f) + (phases(i) - 5 / 9.f) / (1 / 3.f) * vec3(0, 0, 3)));
            target_rots(i) = R * quat(0, vec3(0.f, 1.f, 0.f));
        }
        else
        {
            target_poses(i) = c + R * (length / sqrt(3) * (vec3(- sqrt(3) / 2.f, 0, 3 / 2.f) + (phases(i) - 8 / 9.f) / (1 / 9.f) * vec3(sqrt(3) / 2.f, 0, - 1 / 2.f)));
            target_rots(i) = R * quat(120, vec3(0.f, 1.f, 0.f));
        }
    }
}

struct two_wheeled_car
{
    // 车宽
    float L;
    // 轮子半径
    float R;
    // 车的位置
    vec3 pos;
    // 车的姿态
    quat rot;
    // 左轮速度大小
    float length_vl;
    // 右轮速度大小
    float length_vr;
    float acc_length_vl;
    float acc_length_vr;
    float ek_length_vl;
    float ek_length_vr;

    void vel_control(vec3 target_pos);
    void kinematics(float dt);
};

void two_wheeled_car::vel_control(vec3 target_pos)
{
    // 使用当前位置点与期望位置点的欧氏距离作为车轮平均速度的期望
    float length_vel = length(vec_to_vel(pos, target_pos));
    // 由当前速度方向与（当前位置点指向期望位置点的）方向向量计算出角速度大小
    float length_avel = quat_to_avel(quat(1.f, 0.f, 0.f, 0.f), vec_to_quat(rot * vec3(0.f, 0.f, 1.f), target_pos - pos)).y;
    // float * quat * vec3!!!
    // 计算左右轮期望速度大小
    float target_length_vl = (length_vel - length_avel * L / 2.f);
    float target_length_vr = (length_vel + length_avel * L / 2.f);

    ek_length_vl = ek_length_vl + target_length_vl - length_vl;
    ek_length_vr = ek_length_vr + target_length_vr - length_vr;

    acc_length_vl = 100.f * (1.f * (target_length_vl - length_vl) + 0.05f * ek_length_vl);
    acc_length_vr = 100.f * (1.f * (target_length_vr - length_vr) + 0.05f * ek_length_vr);
}

void two_wheeled_car::kinematics(float dt)
{
    length_vl = length_vl + acc_length_vl * dt;
    length_vr = length_vr + acc_length_vr * dt;
    // 计算实际速度方向
    vec3 vl = length_vl * (rot * vec3(0.f, 0.f, 1.f));
    vec3 vr = length_vr * (rot * vec3(0.f, 0.f, 1.f));
    // 平均速度更新位置
    pos = pos + 0.5f * (vl + vr) * dt;
    // 速度之差得到角速度，更新姿态
    rot = avel_to_quat(cross(rot * vec3(-1, 0.f, 0.f), (vr - vl)) / L / L, dt) * rot;
}

#endif
