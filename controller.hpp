#ifndef ENGINE_CONTROLLER
#define ENGINE_CONTROLLER

#include "vec.hpp"

enum Bait
{
    stand = 0,
    walk  = 1,
    run   = 2,
    jump  = 3,
    squat = 4,
    roll  = 5
};

struct Controller
{
    vec3 *pos;
    vec3 vel;
    vec3 acc;

    vec3 *dir;
    vec2 ang;
    vec2 ang_vel;
    vec2 ang_acc;

    Bait bait;

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

void Controller::pos_pid_control(
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

void Controller::dir_pid_control(
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

void Controller::dir_gamepad_control(vec3 input, float k)
{
    ang.x = circulate_float(ang.x - k * input.x, - PI, PI);
    ang.y = clampf(ang.y + k * input.z, PI / 2, PI - 0.1f);
    *dir = sph_to_dir(ang);
}

void bind_controller(Controller &c, vec3 *pos, vec3 *dir)
{
    c.pos = pos;
    c.dir = dir;
    c.ang = dir_to_sph(*dir);
}

#endif