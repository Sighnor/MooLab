#ifndef ENGINE_MOTION
#define ENGINE_MOTION

#include "array.hpp"
#include "quat.hpp"
#include "spring.hpp"

enum Motion_Bones
{
    RootJoint          = 0,
    lHip          = 1,
    lKnee     = 2,
    lAnkle       = 3,
    lToeJoint      = 4,
    lToeJoint_end       = 5,
    rHip    = 6,
    rKnee      = 7,
    rAnkle     = 8,
    rToeJoint      = 9,
    rToeJoint_end         = 10,
    pelvis_lowerback        = 11,
    lowerback_torso        = 12,
    torso_head          = 13,
    torso_head_end          = 14,
    lTorso_Clavicle  = 15,
    lShoulder       = 16,
    lElbow   = 17,
    lWrist      = 18,
    lWrist_end = 19,
    rTorso_Clavicle      = 20,
    rShoulder  = 21,
    rElbow     = 22,
    rWrist     = 23,
    rWrist_end       = 24,
};

struct BVH_Motion
{
    array1d<int> bone_ids;
    array2d<vec3> bone_positions;
    array2d<quat> bone_rotations;
    array1d<int> bone_parents;

    int nframes() const { return bone_positions.rows; }
    int nbones() const { return bone_positions.cols; }

    // BVH_Motion& operator = (const BVH_Motion& motion) 
    // {
    //     bone_ids = motion.bone_ids;
    //     bone_positions = motion.bone_positions;
    //     bone_rotations = motion.bone_rotations;
    //     bone_parents = motion.bone_parents;

    //     return *this;
    // } 

    BVH_Motion() {}
    BVH_Motion(int nframes, int nbones)
    {
        bone_ids.resize(nbones);
        bone_positions.resize(nframes, nbones);
        bone_rotations.resize(nframes, nbones);
        bone_parents.resize(nbones);

        bone_ids.zero();
        bone_positions.zero();
        bone_rotations.zero();
        bone_parents.zero();
    }
};

void Motion_load(BVH_Motion &motion, const char* file_name)
{
    FILE* f = fopen(file_name, "rb");
    assert(f != NULL);

    array1d_read(motion.bone_ids, f);
    array2d_read(motion.bone_positions, f);
    array2d_read(motion.bone_rotations, f);
    array1d_read(motion.bone_parents, f);

    fclose(f);
}

void Motion_save(BVH_Motion &motion, const char* file_name)
{
    FILE* f = fopen(file_name, "wb");
    assert(f != NULL);

    array1d_write(motion.bone_ids, f);
    array2d_write(motion.bone_positions, f);
    array2d_write(motion.bone_rotations, f);
    array1d_write(motion.bone_parents, f);

    fclose(f);
}

BVH_Motion motion_sub_sequence(const BVH_Motion &motion, int begin, int end)
{
    BVH_Motion res;

    res.bone_ids = motion.bone_ids;
    res.bone_parents = motion.bone_parents;
    res.bone_positions = array2d__sub_sequence(motion.bone_positions, begin, end);
    res.bone_rotations = array2d__sub_sequence(motion.bone_rotations, begin, end);

    return res;
}

BVH_Motion motion_concatenate(const BVH_Motion &motion1, const BVH_Motion &motion2)
{
    BVH_Motion res;

    res.bone_ids = motion1.bone_ids;
    res.bone_parents = motion1.bone_parents;
    res.bone_positions = array2d__concatenate(motion1.bone_positions, motion2.bone_positions);
    res.bone_rotations = array2d__concatenate(motion1.bone_rotations, motion2.bone_rotations);

    return res;
}

void batch_forward_kinematics(
    const BVH_Motion &motion,
    int frame_num, 
    slice1d<vec3> bone_anim_positions, 
    slice1d<quat> bone_anim_rotations)
{
    bone_anim_positions.zero();
    bone_anim_rotations.zero();
    
    for(int i = 0; i < motion.nbones(); i++)
    {
        // Assumes bones are always sorted from root onwards
        int parent_id = motion.bone_parents(i);

        if(parent_id == -1)
        {
            bone_anim_rotations(i) = motion.bone_rotations(frame_num, i);
            bone_anim_positions(i) = motion.bone_positions(frame_num, i);
        }
        else
        {
            bone_anim_rotations(i) = 
                bone_anim_rotations(parent_id) * motion.bone_rotations(frame_num, i);
            bone_anim_positions(i) = 
                bone_anim_positions(parent_id) + 
                bone_anim_rotations(parent_id) * motion.bone_positions(frame_num, i);
        }
    }
}

void decompose_rotation_with_yaxis(
    quat &Ry, 
    quat &Rxz, 
    const quat &rotation)
{
    // vec3 yxz = quat_to_euler_YXZ(rotation);
    // Ry = euler_YXZ_to_quat(yxz.x, 0.f, 0.f);
    // Rxz = euler_YXZ_to_quat(0.f, yxz.y, yxz.z);
    vec2 y_xz = dir_to_angle(rotation * vec3(0.f, 0.f, 1.f));
    Ry = euler_YXZ_to_quat(rad_to_deg(y_xz.x), 0.f, 0.f);
    Rxz = euler_YXZ_to_quat(0.f, rad_to_deg(y_xz.y - PI / 2.f), 0.f);
}

BVH_Motion translation_and_rotation(
    const BVH_Motion &motion, 
    int frame_num, 
    vec3 target_translation_xz, 
    vec3 target_facing_direction_xz)
{
    BVH_Motion res;

    res.bone_ids = motion.bone_ids;
    res.bone_parents = motion.bone_parents;
    res.bone_positions = motion.bone_positions;
    res.bone_rotations = motion.bone_rotations;
  
    vec3 offset = target_translation_xz - motion.bone_positions(frame_num, 0);

    for(int t = 0; t < res.nframes(); t++)
    {
        res.bone_positions(t, 0) = res.bone_positions(t, 0) + offset;
    }

    quat Ry, Rxz;
    decompose_rotation_with_yaxis(Ry, Rxz, res.bone_rotations(frame_num, 0));
    quat invRyxz = inv_quat(Ry * Rxz);

    if(target_facing_direction_xz.x == 0.f && target_facing_direction_xz.z > 0.f)
    {
        Ry = rot_vec_to_quat(vec3(0.f, 1.f, 0.f) * 0.f);
    }
    else if(target_facing_direction_xz.x == 0.f && target_facing_direction_xz.z < 0.f)
    {
        Ry = rot_vec_to_quat(vec3(0.f, 1.f, 0.f) * PI);
    }
    else
    {
        vec3 vector1 = normalize(vec3(0.f, 0.f, 1.f));
        vec3 vector2 = normalize(vec3(target_facing_direction_xz.x, 0.f, target_facing_direction_xz.z));
        vec3 rot_axis = normalize(cross(vector1, vector2));
        float rot_angle = acos(clampf(dot(vector1, vector2), -1.f, 1.f));
        Ry = rot_vec_to_quat(rot_axis * rot_angle);
    }

    quat Q = (Ry * Rxz * invRyxz);

    for(int t = 0; t < res.nframes(); t++)
    {
        res.bone_rotations(t, 0) = Q * res.bone_rotations(t, 0);
        res.bone_positions(t, 0) = res.bone_positions(frame_num, 0) + 
            Q * (res.bone_positions(t, 0) - res.bone_positions(frame_num, 0));
    }

    return res;
}

BVH_Motion blend_two_motions(
    const BVH_Motion &motion1, 
    const BVH_Motion &motion2, 
    array1d<float> alpha)
{
    BVH_Motion res(alpha.size, motion1.nbones());

    for(int i = 0; i < alpha.size; i++)
    {
        float t = float(i) / alpha.size;
        int j = int(t * motion1.nframes());
        int k = int(t * motion2.nframes());

        for(int id = 0; id < res.nbones(); id++)
        {
            res.bone_positions(i, id) = 
                lerp(motion1.bone_positions(j, id), motion2.bone_positions(k, id), alpha(i));

            res.bone_rotations(i, id) = 
                slerp(motion1.bone_rotations(j, id), motion2.bone_rotations(k, id), alpha(i));
        }
    }
    return res; 
}

BVH_Motion build_loop_motion(
    const BVH_Motion &motion, 
    float half_life, 
    float fps)
{
    BVH_Motion res;

    res.bone_ids = motion.bone_ids;
    res.bone_parents = motion.bone_parents;
    res.bone_positions = motion.bone_positions;
    res.bone_rotations = motion.bone_rotations;
    
    for(int ib = 0; ib < res.nbones(); ib++)
    {
        quat avel_begin = 1 / 60.f * res.bone_rotations(0, ib);
        quat avel_end = 1 / 60.f * res.bone_rotations(res.nframes() - 1, ib);
    
        vec3 rot_diff = quat_to_avel(avel_end * inv_quat(avel_begin));
        vec3 avel_diff = quat_to_avel(avel_end) - quat_to_avel(avel_begin);

        vec3 vel1 = res.bone_positions(res.nframes() - 1, ib) - res.bone_positions(res.nframes() - 2, ib);
        vec3 vel2 = res.bone_positions(1, ib) - res.bone_positions(0, ib);

        vec3 pos_diff = vec3(00.f, res.bone_positions(res.nframes() - 1, ib).y - res.bone_positions(0, ib).y, 0.f);
        vec3 vel_diff = (vel1 - vel2) / 60.f;
        
        //旋转差均匀分布到每一帧
        for(int t = 0; t < res.nframes(); t++)
        {
            vec3 offset_rot1 = decay_spring_implicit_damping_rot(
                0.5f * rot_diff, 0.5f * avel_diff, half_life, t / fps);
            vec3 offset_rot2 = decay_spring_implicit_damping_rot(
                -0.5 * rot_diff, -0.5 * avel_diff, half_life, (res.nframes() - t - 1) / fps);
            quat offset_rot = rot_vec_to_quat(offset_rot1 + offset_rot2);

            vec3 offset_pos1 = decay_spring_implicit_damping_pos(
                0.5 * pos_diff, 0.5 * vel_diff, half_life, t / fps);
            vec3 offset_pos2 = decay_spring_implicit_damping_pos(
                -0.5 * pos_diff, -0.5 * vel_diff, half_life, (res.nframes() - t - 1) / fps);
            vec3 offset_pos = offset_pos1 + offset_pos2;

            res.bone_rotations(t, ib) = offset_rot * motion.bone_rotations(t, ib);

            res.bone_positions(t, ib) =  offset_pos + res.bone_positions(t, ib);
        }
    }

    return res;
}

BVH_Motion concatenate_two_motions(
    const BVH_Motion &motion1, 
    const BVH_Motion &motion2,
    int mix_frame,
    int mix_time)
{
    BVH_Motion res1;

    res1.bone_ids = motion1.bone_ids;
    res1.bone_parents = motion1.bone_parents;
    res1.bone_positions = motion1.bone_positions;
    res1.bone_rotations = motion1.bone_rotations;
    
    vec3 pos = res1.bone_positions(mix_frame, 0);
    quat rot = res1.bone_rotations(mix_frame, 0);
    vec3 facing_axis = rot * vec3(0, 0, 1);
    BVH_Motion res2 = translation_and_rotation(motion2, 0, pos, facing_axis);

    BVH_Motion mix_motion0 = motion_sub_sequence(res1, 0, mix_frame);
    BVH_Motion mix_motion1 = motion_sub_sequence(res1, mix_frame, mix_frame + mix_time);
    BVH_Motion mix_motion2 = motion_sub_sequence(res2, 0, mix_time);
    BVH_Motion mix_motion3 = motion_sub_sequence(res2, mix_time, motion2.nframes() - 1);

    array1d<float> alpha(mix_time);
    for(int t = 0; t < mix_time; t++)
    {
        alpha(t) = t / float(mix_time);
    }
    BVH_Motion mix_motion = blend_two_motions(mix_motion1, mix_motion2, alpha);

    BVH_Motion res = motion_concatenate(mix_motion0, mix_motion);
    res = motion_concatenate(res, mix_motion3);
    
    return res;
}

#endif