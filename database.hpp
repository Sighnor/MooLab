#ifndef ENGINE_DATABASE
#define ENGINE_DATABASE

#include "motion.hpp"

struct Database
{
    array1d<int> bone_ids;
    array2d<vec3> bone_positions;
    array2d<vec3> bone_velocities;
    array2d<quat> bone_rotations;
    array2d<vec3> bone_angular_velocities;
    array1d<int> bone_parents;

    int nframes() const { return bone_positions.rows; }
    int nbones() const { return bone_positions.cols; }
};

void Database_load(Database &db, const char* file_name)
{
    FILE* f = fopen(file_name, "rb");
    assert(f != NULL);

    array1d_read(db.bone_ids, f);
    array2d_read(db.bone_positions, f);
    array2d_read(db.bone_velocities, f);
    array2d_read(db.bone_rotations, f);
    array2d_read(db.bone_angular_velocities, f);
    array1d_read(db.bone_parents, f);

    fclose(f);
}

void Database_save(Database &db, const char* file_name)
{
    FILE* f = fopen(file_name, "wb");
    assert(f != NULL);

    array1d_write(db.bone_ids, f);
    array2d_write(db.bone_positions, f);
    array2d_write(db.bone_velocities, f);
    array2d_write(db.bone_rotations, f);
    array2d_write(db.bone_angular_velocities, f);
    array1d_write(db.bone_parents, f);

    fclose(f);
}

void batch_forward_kinematics(
    const Database &db,
    int frame_num,
    slice1d<vec3> bone_anim_positions,
    slice1d<quat> bone_anim_rotations)
{
    bone_anim_positions.zero();
    bone_anim_rotations.zero();

    for (int i = 0; i < db.bone_parents.size; i++)
    {
        // Assumes bones are always sorted from root onwards
        int parent_id = db.bone_parents(i);

        assert(parent_id < i);
        
        if (parent_id == -1)
        {
            bone_anim_rotations(i) = db.bone_rotations(frame_num, i);
            bone_anim_positions(i) = db.bone_positions(frame_num, i);
        }
        else
        {
            bone_anim_rotations(i) = 
                bone_anim_rotations(parent_id) * db.bone_rotations(frame_num, i);
            bone_anim_positions(i) = 
                bone_anim_positions(parent_id) + 
                bone_anim_rotations(parent_id) * db.bone_positions(frame_num, i);
        }
    }
}

#endif