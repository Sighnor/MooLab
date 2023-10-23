#ifndef MOOLAB_SIMULATOR
#define MOOLAB_SIMULATOR

#include "motion.hpp"

struct Simulator
{
    array1d<int> bone_ids;
    array1d<vec3> bone_local_positions;
    array1d<quat> bone_local_rotations;
    array1d<vec3> bone_local_vels;
    array1d<vec3> bone_local_avels;
    array1d<float> bone_masses;
    array1d<float> bone_inertias;
    array1d<float> bone_local_forces;
    array1d<vec3> bone_local_torgues;
    array1d<int> bone_parents;
    array1d<vec3> bone_anim_positions;
    array1d<quat> bone_anim_rotations;

    int nbones() const { return bone_local_positions.size; }

    void bind_simulator(const BVH_Motion &motion);
    void simulate_gravity();
    void simulate(float dt);
    void batch_forward_kinematics_full();
};

void Simulator::bind_simulator(const BVH_Motion &motion)
{
    int size = motion.nbones();
    bone_ids = motion.bone_ids;
    bone_local_positions = motion.bone_local_positions(0);
    bone_local_rotations = motion.bone_local_rotations(0);
    bone_local_vels.resize(size);
    bone_local_avels.resize(size);
    bone_local_avels.zero();
    bone_masses.resize(size);
    bone_inertias.resize(size);
    bone_local_forces.resize(size);
    bone_local_torgues.resize(size);
    bone_local_torgues.zero();
    bone_parents = motion.bone_parents;
    bone_anim_positions.resize(size);
    bone_anim_rotations.resize(size);
}

void Simulator::simulate_gravity()
{
    for(int i = 22; i < nbones(); i++)
    {
        int parent_id = bone_parents(i);

        if(((bone_anim_positions(parent_id) + bone_anim_positions(i)) / 2).y > 0.f)
        {
            vec3 r = 0.5f * (bone_anim_rotations(parent_id) * bone_local_positions(i));
            vec3 g = vec3(0.f, 1.f, 0.f);
            bone_local_torgues(parent_id) = bone_anim_rotations(parent_id) * cross(r, g);
        }
    }
}

void Simulator::simulate(float dt = 0.0166667f)
{
    for(int i = 1; i < nbones(); i++)
    {
        bone_local_avels(i) = 5.f * bone_local_torgues(i);
        bone_local_rotations(i) = avel_to_quat(bone_local_avels(i), dt) * bone_local_rotations(i);
    }
}

void Simulator::batch_forward_kinematics_full()
{
    bone_anim_positions.zero();
    bone_anim_rotations.zero();
    
    for(int i = 0; i < nbones(); i++)
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