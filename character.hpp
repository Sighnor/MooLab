#ifndef ENGINE_CHARACTER
#define ENGINE_CHARACTER

#include "array.hpp"
#include "mesh.hpp"
#include "quat.hpp"

enum Bones
{
    Bone_Entity        = 0,
    Bone_Hips          = 1,
    Bone_LeftUpLeg     = 2,
    Bone_LeftLeg       = 3,
    Bone_LeftFoot      = 4,
    Bone_LeftToe       = 5,
    Bone_RightUpLeg    = 6,
    Bone_RightLeg      = 7,
    Bone_RightFoot     = 8,
    Bone_RightToe      = 9,
    Bone_Spine         = 10,
    Bone_Spine1        = 11,
    Bone_Spine2        = 12,
    Bone_Neck          = 13,
    Bone_Head          = 14,
    Bone_LeftShoulder  = 15,
    Bone_LeftArm       = 16,
    Bone_LeftForeArm   = 17,
    Bone_LeftHand      = 18,
    Bone_RightShoulder = 19,
    Bone_RightArm      = 20,
    Bone_RightForeArm  = 21,
    Bone_RightHand     = 22,
};

struct Character
{
    array1d<vec3> all_rest_positions;
    array1d<vec3> all_rest_normals;
    array1d<vec2> texcoords;
    array1d<unsigned short> indices;

    array2d<float> bone_weights;
    array2d<unsigned short> bone_weights_ids;

    array1d<vec3> bone_rest_positions;
    array1d<quat> bone_rest_rotations;
};

void character_load(Character& character, const char* filename)
{
    FILE* f = fopen(filename, "rb");
    assert(f != NULL);
    
    array1d_read(character.all_rest_positions, f);
    array1d_read(character.all_rest_normals, f);
    array1d_read(character.texcoords, f);
    array1d_read(character.indices, f);
    
    array2d_read(character.bone_weights, f);
    array2d_read(character.bone_weights_ids, f);
    
    array1d_read(character.bone_rest_positions, f);
    array1d_read(character.bone_rest_rotations, f);
    
    fclose(f);
}

void linear_blend_positions(
    slice1d<vec3> all_anim_positions,
    const slice1d<vec3> bone_anim_positions,
    const slice1d<quat> bone_anim_rotations,
    const slice1d<vec3> all_rest_positions,
    const slice1d<vec3> bone_rest_positions,
    const slice1d<quat> bone_rest_rotations,
    const slice2d<unsigned short> bone_weights_ids,
    const slice2d<float> bone_weights)
{
    all_anim_positions.zero();

    for(int i = 0; i < all_anim_positions.size; i++)
    {
        for(int j = 0; j < bone_weights_ids.cols; j++)
        {
            if(bone_weights(i, j) > 0.f)
            {
                int id = bone_weights_ids(i, j);

                vec3 position = bone_anim_rotations(id) * (inv_quat(bone_rest_rotations(id)) * (all_rest_positions(i) - bone_rest_positions(id)))
                                + bone_anim_positions(id);

                all_anim_positions(i) = all_anim_positions(i) + bone_weights(i, j) * position;
            }
        }
    }
}

void linear_blend_normals(
    slice1d<vec3> all_anim_normals,
    const slice1d<quat> bone_anim_rotations,
    const slice1d<vec3> all_rest_normals,
    const slice1d<quat> bone_rest_rotations,
    const slice2d<unsigned short> bone_weights_ids,
    const slice2d<float> bone_weights)
{
    all_anim_normals.zero();

    for(int i = 0; i < all_anim_normals.size; i++)
    {
        for(int j = 0; j < bone_weights_ids.cols; j++)
        {
            if(bone_weights(i, j) > 0.f)
            {
                int id = bone_weights_ids(i, j);

                vec3 normal = bone_anim_rotations(id) * (inv_quat(bone_rest_rotations(id)) * all_rest_normals(i));

                all_anim_normals(i) = all_anim_normals(i) + bone_weights(i, j) * normal;
            }
        }
    }
}

MooMesh make_character_rest_mesh(const Character &character)
{
    MooMesh mesh;

    mesh.vertex_count = character.all_rest_positions.size;
    mesh.triangle_count = character.indices.size / 3;
    mesh.vertices = (float *)malloc(character.all_rest_positions.size * 3 * sizeof(float));
    mesh.normals = (float *)malloc(character.all_rest_normals.size * 3 * sizeof(float));
    mesh.texcoords = (float *)malloc(character.texcoords.size * 2 * sizeof(float));
    mesh.indices = (unsigned short *)malloc(character.indices.size * sizeof(unsigned short));

    memcpy(mesh.vertices, character.all_rest_positions.data, character.all_rest_positions.size * 3 * sizeof(float));
    memcpy(mesh.normals, character.all_rest_normals.data, character.all_rest_normals.size * 3 * sizeof(float));
    memcpy(mesh.texcoords, character.texcoords.data, character.texcoords.size * 2 * sizeof(float));
    memcpy(mesh.indices, character.indices.data, character.indices.size * sizeof(unsigned short));

    return mesh;
}

void deform_character_anim_mesh(
    const Character &character, 
    const slice1d<vec3> bone_anim_positions,
    const slice1d<quat> bone_anim_rotations,
    MooMesh &mesh)
{
    linear_blend_positions(
        slice1d<vec3>(mesh.vertex_count, (vec3 *)mesh.vertices),
        bone_anim_positions,
        bone_anim_rotations,
        character.all_rest_positions,
        character.bone_rest_positions,
        character.bone_rest_rotations,
        character.bone_weights_ids,
        character.bone_weights);
    linear_blend_normals(
        slice1d<vec3>(mesh.vertex_count, (vec3 *)mesh.normals),
        bone_anim_rotations,
        character.all_rest_normals,
        character.bone_rest_rotations,
        character.bone_weights_ids,
        character.bone_weights);
}

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

    float lab = length(ab);
    float lbc = length(bc);
    float lac = length(ac);
    float lat = length(at);

    vec3 axis0 = normalize(cross(ac, ab));
    vec3 axis1 = normalize(cross(bc, ab));
    vec3 axis2 = normalize(cross(ac, at));

    float thetaa0 = acos(clampf((lab * lab + lac * lac - lbc * lbc) / (2.f * lab * lac), -1.f, 1.f));
    float thetaa1 = acos(clampf((lab * lab + lat * lat - lbc * lbc) / (2.f * lab * lat), -1.f, 1.f));
    float thetab0 = acos(clampf((lab * lab + lbc * lbc - lac * lac) / (2.f * lab * lbc), -1.f, 1.f));
    float thetab1 = acos(clampf((lab * lab + lbc * lbc - lat * lat) / (2.f * lab * lbc), -1.f, 1.f));
    float thetat = acos(dot(ac, at) / (lac * lat));

    quat R0 = quat(rad_to_deg(thetaa1 - thetaa0), axis0);
    quat R1 = quat(rad_to_deg(thetab1 - thetab0), axis1);
    quat R2 = quat(rad_to_deg(thetat), axis2);

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

    float lab = length(ab);
    float lbc = length(bc);
    float lac = length(ac);
    float lat = length(at);

    vec3 axis0 = normalize(cross(ac, ab));
    vec3 axis1 = normalize(cross(bc, ab));
    vec3 axis2 = normalize(cross(ac, at));

    float thetaa0 = acos(clampf((lab * lab + lac * lac - lbc * lbc) / (2.f * lab * lac), -1.f, 1.f));
    float thetaa1 = acos(clampf((lab * lab + lat * lat - lbc * lbc) / (2.f * lab * lat), -1.f, 1.f));
    float thetab0 = acos(clampf((lab * lab + lbc * lbc - lac * lac) / (2.f * lab * lbc), -1.f, 1.f));
    float thetab1 = acos(clampf((lab * lab + lbc * lbc - lat * lat) / (2.f * lab * lbc), -1.f, 1.f));
    float thetat = acos(dot(ac, at) / (lac * lat));

    quat R0 = quat(rad_to_deg(thetaa1 - thetaa0), axis0);
    quat R1 = quat(rad_to_deg(thetab1 - thetab0), axis1);
    quat R2 = quat(rad_to_deg(thetat), axis2);

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

#endif