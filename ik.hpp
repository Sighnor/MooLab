#ifndef ENGINE_IK
#define ENGINE_IK

#include "array.hpp"
#include "quat.hpp"

void IK_two_bones(
    vec3 &bone0_anim_position, 
    vec3 &bone1_anim_position, 
    vec3 &bone2_anim_position,
    quat &bone0_anim_rotation, 
    quat &bone1_anim_rotation, 
    quat &bone2_anim_rotation, 
    vec3 target_position)
{
    vec3 ab = bone1_anim_position - bone0_anim_position;
    vec3 bc = bone2_anim_position - bone1_anim_position;
    vec3 ac = bone2_anim_position - bone0_anim_position;
    vec3 at = target_position - bone0_anim_position;

    float l_ab = length(ab);
    float l_bc = length(bc);
    float l_ac = length(ac);
    float l_at = length(at);

    vec3 axis0 = normalize(cross(ac, ab));
    vec3 axis1 = normalize(cross(bc, ab));
    vec3 axis2 = normalize(cross(ac, at));

    float theta_a0 = acos(clampf((l_ab * l_ab + l_ac * l_ac - l_bc * l_bc) / (2.f * l_ab * l_ac), -1.f, 1.f));
    float theta_a1 = acos(clampf((l_ab * l_ab + l_at * l_at - l_bc * l_bc) / (2.f * l_ab * l_at), -1.f, 1.f));
    float theta_b0 = acos(clampf((l_ab * l_ab + l_bc * l_bc - l_ac * l_ac) / (2.f * l_ab * l_bc), -1.f, 1.f));
    float theta_b1 = acos(clampf((l_ab * l_ab + l_bc * l_bc - l_at * l_at) / (2.f * l_ab * l_bc), -1.f, 1.f));
    float theta_t = acos(clampf(dot(ac, at) / (l_ac * l_at), -1.f, 1.f));

    quat R0 = quat(rad_to_deg(theta_a1 - theta_a0), axis0);
    quat R1 = quat(rad_to_deg(theta_b1 - theta_b0), axis1);
    quat R2 = quat(rad_to_deg(theta_t), axis2);

    vec3 offset01 = inv_quat(bone0_anim_rotation) * (bone1_anim_position - bone0_anim_position);
    vec3 offset12 = inv_quat(bone1_anim_rotation) * (bone2_anim_position - bone1_anim_position);
    quat rotation01 = inv_quat(bone0_anim_rotation) * bone1_anim_rotation;
    quat rotation12 = inv_quat(bone1_anim_rotation) * bone2_anim_rotation;

    bone0_anim_position = bone0_anim_position;
    bone0_anim_rotation = R0 * R2 * bone0_anim_rotation;
    bone1_anim_position = bone0_anim_position + bone0_anim_rotation * offset01;
    bone1_anim_rotation = R1 * bone0_anim_rotation * rotation01;
    bone2_anim_position = bone1_anim_position + bone1_anim_rotation * offset12;
    bone2_anim_rotation = bone1_anim_rotation * rotation12;
}

void IK_two_bones(
    vec3 &bone0_anim_position, 
    vec3 &bone1_anim_position, 
    vec3 &bone2_anim_position,
    vec3 &bone3_anim_position,
    quat &bone0_anim_rotation, 
    quat &bone1_anim_rotation, 
    quat &bone2_anim_rotation, 
    quat &bone3_anim_rotation, 
    vec3 target_position)
{
    vec3 ab = bone1_anim_position - bone0_anim_position;
    vec3 bc = bone2_anim_position - bone1_anim_position;
    vec3 ac = bone2_anim_position - bone0_anim_position;
    vec3 at = target_position - bone0_anim_position;

    float l_ab = length(ab);
    float l_bc = length(bc);
    float l_ac = length(ac);
    float l_at = length(at);

    vec3 axis0 = normalize(cross(ac, ab));
    vec3 axis1 = normalize(cross(bc, ab));
    vec3 axis2 = normalize(cross(ac, at));

    float theta_a0 = acos(clampf((l_ab * l_ab + l_ac * l_ac - l_bc * l_bc) / (2.f * l_ab * l_ac), -1.f, 1.f));
    float theta_a1 = acos(clampf((l_ab * l_ab + l_at * l_at - l_bc * l_bc) / (2.f * l_ab * l_at), -1.f, 1.f));
    float theta_b0 = acos(clampf((l_ab * l_ab + l_bc * l_bc - l_ac * l_ac) / (2.f * l_ab * l_bc), -1.f, 1.f));
    float theta_b1 = acos(clampf((l_ab * l_ab + l_bc * l_bc - l_at * l_at) / (2.f * l_ab * l_bc), -1.f, 1.f));
    float theta_t = acos(clampf(dot(ac, at) / (l_ac * l_at), -1.f, 1.f));

    quat R0 = quat(rad_to_deg(theta_a1 - theta_a0), axis0);
    quat R1 = quat(rad_to_deg(theta_b1 - theta_b0), axis1);
    quat R2 = quat(rad_to_deg(theta_t), axis2);

    vec3 offset01 = inv_quat(bone0_anim_rotation) * (bone1_anim_position - bone0_anim_position);
    vec3 offset12 = inv_quat(bone1_anim_rotation) * (bone2_anim_position - bone1_anim_position);
    vec3 offset23 = inv_quat(bone2_anim_rotation) * (bone3_anim_position - bone2_anim_position);
    quat rotation01 = inv_quat(bone0_anim_rotation) * bone1_anim_rotation;
    quat rotation12 = inv_quat(bone1_anim_rotation) * bone2_anim_rotation;
    quat rotation23 = inv_quat(bone2_anim_rotation) * bone3_anim_rotation;

    bone0_anim_position = bone0_anim_position;
    bone0_anim_rotation = R0 * R2 * bone0_anim_rotation;
    bone1_anim_position = bone0_anim_position + bone0_anim_rotation * offset01;
    bone1_anim_rotation = R1 * bone0_anim_rotation * rotation01;
    bone2_anim_position = bone1_anim_position + bone1_anim_rotation * offset12;
    bone2_anim_rotation = bone1_anim_rotation * rotation12;
    bone3_anim_position = bone2_anim_position + bone2_anim_rotation * offset23;
    bone3_anim_rotation = bone2_anim_rotation * rotation23;
}

void FABRIK_one_end(
    const slice1d<vec3> bone_local_positions, 
    slice1d<quat> bone_local_rotations, 
    vec3 target0, 
    int id0, 
    float eps = 0.0001f, 
    int max = 10)
{
    int size = bone_local_positions.size;
    int t = 0;
    vec3 ori = bone_local_positions(0);

    array1d<vec3> bone_anim_positions(size);
    array1d<quat> bone_anim_rotations(size);
    array1d<float> lengths(size);
    
    bone_anim_positions(0) = ori;
    lengths(0) = 0.f;
    for(int i = 1; i < size; i++)
    {
        bone_anim_positions(i) = bone_anim_positions(i - 1) + bone_local_positions(i);
        lengths(i) = length(bone_local_positions(i));
    }

    while (length(bone_anim_positions(id0) - target0) > eps && t < max)
    {
        //Backward
        bone_anim_positions(id0) = target0;
        for(int i = id0 - 1; i > -1; i--)
        {
            vec3 v0 = normalize(bone_anim_positions(i) - bone_anim_positions(i + 1));
            bone_anim_positions(i) = bone_anim_positions(i + 1) + lengths(i + 1) * v0;
        }
        //Forward
        bone_anim_positions(0) = ori;
        for(int i = 1; i < id0 + 1; i++)
        {
            vec3 v0 = normalize(bone_anim_positions(i) - bone_anim_positions(i - 1));
            bone_anim_positions(i) = bone_anim_positions(i - 1) + lengths(i) * v0;
        }
        t++;
    }

    for(int i = 1; i < id0 + 1; i++)
    {
        vec3 v0 = normalize(bone_local_positions(i));
        vec3 v1 = normalize(bone_anim_positions(i) - bone_anim_positions(i - 1));
        vec3 axis = normalize(cross(v0, v1));
        float phi = rad_to_deg(acos(clampf(dot(v0, v1), -1.f, 1.f)));
        bone_anim_rotations(i - 1) = quat(phi, axis);
    }

    bone_local_rotations(0) = bone_anim_rotations(0);
    for(int i = 1; i < id0; i++)
    {
        bone_local_rotations(i) = inv_quat(bone_anim_rotations(i - 1)) * bone_anim_rotations(i);
    }
}

void FABRIK_one_end_constraints(
    const slice1d<vec3> bone_local_positions, 
    slice1d<quat> bone_local_rotations, 
    vec3 target0, 
    int id0, 
    float constraint_angle, 
    float eps = 0.0001f, 
    int max = 10)
{
    int size = bone_local_positions.size;
    int t = 0;
    vec3 ori = bone_local_positions(0);

    array1d<vec3> bone_anim_positions(size);
    array1d<quat> bone_anim_rotations(size);
    array1d<float> lengths(size);
    
    bone_anim_positions(0) = ori;
    lengths(0) = 0.f;
    for(int i = 1; i < size; i++)
    {
        bone_anim_positions(i) = bone_anim_positions(i - 1) + bone_local_positions(i);
        lengths(i) = length(bone_local_positions(i));
    }

    while (length(bone_anim_positions(id0) - target0) > eps && t < max)
    {
        //Backward
        bone_anim_positions(id0) = target0;
        vec3 v0 = normalize(bone_anim_positions(id0 - 1) - bone_anim_positions(id0));
        bone_anim_positions(id0 - 1) = bone_anim_positions(id0) + lengths(id0) * v0;
        for(int i = id0 - 2; i > -1; i--)
        {
            vec3 v1 = normalize(bone_anim_positions(i + 1) - bone_anim_positions(i + 2));
            vec3 v2 = normalize(bone_anim_positions(i) - bone_anim_positions(i + 1));
            if(acos(clampf(dot(v1, v2), -1.f, 1.f)) > deg_to_rad(constraint_angle))
            {
                vec3 axis = normalize(cross(v1, v2));
                vec3 v3 = quat(constraint_angle, axis) * v1;
                bone_anim_positions(i) = bone_anim_positions(i + 1) + lengths(i + 1) * v3;
            }
            else
            {
                vec3 v3 = normalize(bone_anim_positions(i) - bone_anim_positions(i + 1));
                bone_anim_positions(i) = bone_anim_positions(i + 1) + lengths(i + 1) * v3;
            }
        }
        //Forward
        bone_anim_positions(0) = ori;
        vec3 v4 = normalize(bone_anim_positions(1) - bone_anim_positions(0));
        bone_anim_positions(1) = bone_anim_positions(0) + lengths(1) * v4;
        for(int i = 2; i < id0 + 1; i++)
        {
            vec3 v1 = normalize(bone_anim_positions(i - 1) - bone_anim_positions(i - 2));
            vec3 v2 = normalize(bone_anim_positions(i) - bone_anim_positions(i - 1));
            if(acos(clampf(dot(v1, v2), -1.f, 1.f)) > deg_to_rad(constraint_angle))
            {
                vec3 axis = normalize(cross(v1, v2));
                vec3 v3 = quat(constraint_angle, axis) * v1;
                bone_anim_positions(i) = bone_anim_positions(i - 1) + lengths(i) * v3;
            }
            else
            {
                vec3 v3 = normalize(bone_anim_positions(i) - bone_anim_positions(i - 1));
                bone_anim_positions(i) = bone_anim_positions(i - 1) + lengths(i) * v3;
            }
        }
        t++;
    }

    for(int i = 1; i < id0 + 1; i++)
    {
        vec3 v0 = normalize(bone_local_positions(i));
        vec3 v1 = normalize(bone_anim_positions(i) - bone_anim_positions(i - 1));
        vec3 axis = normalize(cross(v0, v1));
        float phi = rad_to_deg(acos(clampf(dot(v0, v1), -1.f, 1.f)));
        bone_anim_rotations(i - 1) = quat(phi, axis);
    }

    bone_local_rotations(0) = bone_anim_rotations(0);
    for(int i = 1; i < id0; i++)
    {
        bone_local_rotations(i) = inv_quat(bone_anim_rotations(i - 1)) * bone_anim_rotations(i);
    }
}

void FABRIK_two_end(
    const slice1d<vec3> bone_local_positions, 
    slice1d<quat> bone_local_rotations, 
    vec3 target0, 
    vec3 target1, 
    int id0, 
    int id1, 
    float eps = 0.0002f, 
    int max = 20)
{
    int size = bone_local_positions.size;
    int t = 0;
    vec3 ori = bone_local_positions(0);

    array1d<vec3> bone_anim_positions(size);
    array1d<quat> bone_anim_rotations(size);
    array1d<float> lengths(size);
    
    bone_anim_positions(0) = ori;
    lengths(0) = 0.f;
    for(int i = 1; i < size; i++)
    {
        bone_anim_positions(i) = bone_anim_positions(i - 1) + bone_local_positions(i);
        lengths(i) = length(bone_local_positions(i));
    }

    while (length(bone_anim_positions(id0) - target0) + length(bone_anim_positions(id1) - target1) > eps && t < max)
    {
        //Backward
        bone_anim_positions(id0) = target0;
        for(int i = id0 - 1; i > id1 - 1; i--)
        {
            vec3 v0 = normalize(bone_anim_positions(i) - bone_anim_positions(i + 1));
            bone_anim_positions(i) = bone_anim_positions(i + 1) + lengths(i + 1) * v0;
        }
        bone_anim_positions(id1) = (bone_anim_positions(id1) + target1) / 2.f;
        for(int i = id1 - 1; i > -1; i--)
        {
            vec3 v0 = normalize(bone_anim_positions(i) - bone_anim_positions(i + 1));
            bone_anim_positions(i) = bone_anim_positions(i + 1) + lengths(i + 1) * v0;
        }
        //Forward
        bone_anim_positions(0) = ori;
        for(int i = 1; i < id0 + 1; i++)
        {
            vec3 v0 = normalize(bone_anim_positions(i) - bone_anim_positions(i - 1));
            bone_anim_positions(i) = bone_anim_positions(i - 1) + lengths(i) * v0;
        }
        t++;
    }

    for(int i = 1; i < id0 + 1; i++)
    {
        vec3 v0 = normalize(bone_local_positions(i));
        vec3 v1 = normalize(bone_anim_positions(i) - bone_anim_positions(i - 1));
        vec3 axis = normalize(cross(v0, v1));
        float phi = rad_to_deg(acos(clampf(dot(v0, v1), -1.f, 1.f)));
        bone_anim_rotations(i - 1) = quat(phi, axis);
    }

    bone_local_rotations(0) = bone_anim_rotations(0);
    for(int i = 1; i < id0; i++)
    {
        bone_local_rotations(i) = inv_quat(bone_anim_rotations(i - 1)) * bone_anim_rotations(i);
    }
}

void FABRIK_two_end_constraints(
    const slice1d<vec3> bone_local_positions, 
    slice1d<quat> bone_local_rotations, 
    vec3 target0, 
    vec3 target1, 
    int id0, 
    int id1, 
    float constraint_angle, 
    float eps = 0.0002f, 
    int max = 20)
{
    int size = bone_local_positions.size;
    int t = 0;
    vec3 ori = bone_local_positions(0);

    array1d<vec3> bone_anim_positions(size);
    array1d<quat> bone_anim_rotations(size);
    array1d<float> lengths(size);
    
    bone_anim_positions(0) = ori;
    lengths(0) = 0.f;
    for(int i = 1; i < size; i++)
    {
        bone_anim_positions(i) = bone_anim_positions(i - 1) + bone_local_positions(i);
        lengths(i) = length(bone_local_positions(i));
    }

    while (length(bone_anim_positions(id0) - target0) + length(bone_anim_positions(id1) - target1) > eps && t < max)
    {
        //Backward
        bone_anim_positions(id0) = target0;
        vec3 v0 = normalize(bone_anim_positions(id0 - 1) - bone_anim_positions(id0));
        bone_anim_positions(id0 - 1) = bone_anim_positions(id0) + lengths(id0) * v0;
        for(int i = id0 - 2; i > id1 - 1; i--)
        {
            vec3 v1 = normalize(bone_anim_positions(i + 1) - bone_anim_positions(i + 2));
            vec3 v2 = normalize(bone_anim_positions(i) - bone_anim_positions(i + 1));
            if(acos(clampf(dot(v1, v2), -1.f, 1.f)) > deg_to_rad(constraint_angle))
            {
                vec3 axis = normalize(cross(v1, v2));
                vec3 v3 = quat(constraint_angle, axis) * v1;
                bone_anim_positions(i) = bone_anim_positions(i + 1) + lengths(i + 1) * v3;
            }
            else
            {
                vec3 v3 = normalize(bone_anim_positions(i) - bone_anim_positions(i + 1));
                bone_anim_positions(i) = bone_anim_positions(i + 1) + lengths(i + 1) * v3;
            }
        }
        bone_anim_positions(id1) = (bone_anim_positions(id1) + target1) / 2.f;
        for(int i = id1 - 1; i > -1; i--)
        {
            vec3 v1 = normalize(bone_anim_positions(i + 1) - bone_anim_positions(i + 2));
            vec3 v2 = normalize(bone_anim_positions(i) - bone_anim_positions(i + 1));
            if(acos(clampf(dot(v1, v2), -1.f, 1.f)) > deg_to_rad(constraint_angle))
            {
                vec3 axis = normalize(cross(v1, v2));
                vec3 v3 = quat(constraint_angle, axis) * v1;
                bone_anim_positions(i) = bone_anim_positions(i + 1) + lengths(i + 1) * v3;
            }
            else
            {
                vec3 v3 = normalize(bone_anim_positions(i) - bone_anim_positions(i + 1));
                bone_anim_positions(i) = bone_anim_positions(i + 1) + lengths(i + 1) * v3;
            }
        }
        //Forward
        bone_anim_positions(0) = ori;
        vec3 v4 = normalize(bone_anim_positions(1) - bone_anim_positions(0));
        bone_anim_positions(1) = bone_anim_positions(0) + lengths(1) * v4;
        for(int i = 2; i < id0 + 1; i++)
        {
            vec3 v1 = normalize(bone_anim_positions(i - 1) - bone_anim_positions(i - 2));
            vec3 v2 = normalize(bone_anim_positions(i) - bone_anim_positions(i - 1));
            if(acos(clampf(dot(v1, v2), -1.f, 1.f)) > deg_to_rad(constraint_angle))
            {
                vec3 axis = normalize(cross(v1, v2));
                vec3 v3 = quat(constraint_angle, axis) * v1;
                bone_anim_positions(i) = bone_anim_positions(i - 1) + lengths(i) * v3;
            }
            else
            {
                vec3 v3 = normalize(bone_anim_positions(i) - bone_anim_positions(i - 1));
                bone_anim_positions(i) = bone_anim_positions(i - 1) + lengths(i) * v3;
            }
        }
        t++;
    }

    for(int i = 1; i < id0 + 1; i++)
    {
        vec3 v0 = normalize(bone_local_positions(i));
        vec3 v1 = normalize(bone_anim_positions(i) - bone_anim_positions(i - 1));
        vec3 axis = normalize(cross(v0, v1));
        float phi = rad_to_deg(acos(clampf(dot(v0, v1), -1.f, 1.f)));
        bone_anim_rotations(i - 1) = quat(phi, axis);
    }

    bone_local_rotations(0) = bone_anim_rotations(0);
    for(int i = 1; i < id0; i++)
    {
        bone_local_rotations(i) = inv_quat(bone_anim_rotations(i - 1)) * bone_anim_rotations(i);
    }
}

#endif