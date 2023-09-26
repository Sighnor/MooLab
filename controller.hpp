#ifndef MOOLAB_CONTROLLER
#define MOOLAB_CONTROLLER

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

#endif