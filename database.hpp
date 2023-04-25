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
    FILE* f = fopen(file_name, "rb");
    assert(f != NULL);

    array1d_write(db.bone_ids, f);
    array2d_write(db.bone_positions, f);
    array2d_write(db.bone_velocities, f);
    array2d_write(db.bone_rotations, f);
    array2d_write(db.bone_angular_velocities, f);
    array1d_write(db.bone_parents, f);

    fclose(f);
}

void forward_kinematics(
    vec3& bone_position,
    quat& bone_rotation,
    const array1d<vec3> bone_positions,
    const array1d<quat> bone_rotations,
    const array1d<int> bone_parents,
    const int bone_id)
{
    if (bone_parents(bone_id) != -1)
    {
        vec3 parent_position;
        quat parent_rotation;
        
        forward_kinematics(
            parent_position,
            parent_rotation,
            bone_positions,
            bone_rotations,
            bone_parents,
            bone_parents(bone_id));
        
        bone_position = (parent_rotation * bone_positions(bone_id)) + parent_position;
        bone_rotation = (parent_rotation * bone_rotations(bone_id));
    }
    else
    {
        bone_position = bone_positions(bone_id);
        bone_rotation = bone_rotations(bone_id); 
    }
}

#endif