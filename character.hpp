#ifndef ENGINE_CHARACTER
#define ENGINE_CHARACTER

#include "array.hpp"
#include "mesh.hpp"
#include "quat.hpp"

struct Character
{
    array1d<vec3> all_rest_positions;
    array1d<vec3> all_rest_normals;
    array1d<vec2> texcoords;
    array1d<int> indices;

    array2d<int> bone_weights_ids;
    array2d<float> bone_weights;

    array1d<vec3> bone_rest_positions;
    array1d<quat> bone_rest_rotations;
};

void linear_blend_positions(
    slice1d<vec3> all_anim_positions,
    const slice1d<vec3> bone_anim_positions,
    const slice1d<quat> bone_anim_rotations,
    const slice1d<vec3> all_rest_positions,
    const slice1d<vec3> bone_rest_positions,
    const slice1d<quat> bone_rest_rotations,
    const slice2d<int> bone_weights_ids,
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

                vec3 position = bone_anim_rotations(i) * (inv_quat(bone_rest_rotations(i)) * (all_rest_positions(i) - bone_rest_positions(id)))
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
    const slice2d<int> bone_weights_ids,
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

                vec3 position = bone_anim_rotations(i) * (inv_quat(bone_rest_rotations(i)) * all_rest_normals(i));

                all_anim_normals(i) = all_anim_normals(i) + bone_weights(i, j) * position;
            }
        }
    }
}

Mesh make_character_rest_mesh(const Character &character)
{
    Mesh mesh;

    mesh.vertex_count = character.all_rest_positions.size;
    mesh.triangle_count = character.indices.size;
    mesh.vertices = (float *)malloc(character.all_rest_positions.size * 3 * sizeof(float));
    mesh.normals = (float *)malloc(character.all_rest_normals.size * 3 * sizeof(float));
    mesh.texcoords = (float *)malloc(character.texcoords.size * 2 * sizeof(float));
    mesh.indices = (int *)malloc(character.indices.size * sizeof(int));

    memcpy(mesh.vertices, character.all_rest_positions.data, character.all_rest_positions.size * 3 * sizeof(float));
    memcpy(mesh.normals, character.all_rest_normals.data, character.all_rest_normals.size * 3 * sizeof(float));
    memcpy(mesh.texcoords, character.texcoords.data, character.texcoords.size * 2 * sizeof(float));
    memcpy(mesh.indices, character.indices.data, character.indices.size * sizeof(int));

    return mesh;
}

void deform_character_anim_mesh(
    const Character &character, 
    const array1d<vec3> bone_anim_positions,
    const array1d<quat> bone_anim_rotations,
    Mesh &mesh)
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

#endif