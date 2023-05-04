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
    vec3 ang_vel;
    vec3 ang_acc;

    Bait bait;

    void pos_pid_control(
            float p_const, 
            float i_const, 
            float d_const, 
            float k, 
            float dt, 
            vec3 input);
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

void dir_pid_control(Controller &c, float p_const, float i_const, float d_const, float k, float dt, vec3 input)
{
    static vec3 dir_ek_0 = vec3(0.f);
    static vec3 dir_ek_1 = vec3(0.f);
    static vec3 dir_ek_2 = vec3(0.f);

    dir_ek_2 = dir_ek_1;
    dir_ek_1 = dir_ek_0;
    dir_ek_0 = input - c.vel;

    c.acc = (p_const * (dir_ek_0 - dir_ek_1) + i_const * dir_ek_0 - d_const * (dir_ek_0 - 2 * dir_ek_1 + dir_ek_2));
    c.vel = c.vel + c.acc * dt;
    *c.dir = *c.dir + c.vel * dt;
}

void bind_controller(Controller &c, vec3 *pos, vec3 *dir)
{
    c.pos = pos;
    c.dir = dir;
}


#endif