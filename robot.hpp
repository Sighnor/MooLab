#ifndef ENGINE_ROBOT
#define ENGINE_ROBOT

#include "character.hpp"

struct Robot
{
    array1d<int> bone_ids;
    array1d<vec3> bone_offset_positions;
    array1d<vec3> bone_rest_positions;
    array1d<quat> bone_rest_rotations;
    array1d<int> bone_parents;

    array1d<vec3> all_rest_positions;
    array1d<vec3> all_rest_normals;
    array1d<vec2> texcoords;
    array1d<unsigned short> indices;

    array2d<float> bone_weights;
    array2d<unsigned short> bone_weights_ids;

    int nbones() const { return bone_rest_positions.size; }
    void resize_bone(int bones_size)
    {
        bone_ids.resize(bones_size);
        bone_offset_positions.resize(bones_size);
        bone_rest_positions.resize(bones_size);
        bone_rest_rotations.resize(bones_size);
        bone_parents.resize(bones_size);
    }
    void resize_mesh(int triangles_size, int vertices_size) 
    { 
        all_rest_positions.resize(vertices_size);
        all_rest_normals.resize(vertices_size);
        texcoords.resize(vertices_size);
        indices.resize(3 * triangles_size);

        bone_weights.resize(vertices_size, 4);
        bone_weights_ids.resize(vertices_size, 4);
    }
};

void robot_load(Robot& robot, const char* filename)
{
    FILE* f = fopen(filename, "rb");
    assert(f != NULL);

    array1d_read(robot.bone_ids, f);
    array1d_read(robot.bone_offset_positions, f);
    array1d_read(robot.bone_rest_positions, f);
    array1d_read(robot.bone_rest_rotations, f);
    array1d_read(robot.bone_parents, f);
    
    array1d_read(robot.all_rest_positions, f);
    array1d_read(robot.all_rest_normals, f);
    array1d_read(robot.texcoords, f);
    array1d_read(robot.indices, f);
    
    array2d_read(robot.bone_weights, f);
    array2d_read(robot.bone_weights_ids, f);
    
    fclose(f);
}

void robot_save(Robot& robot, const char* filename)
{
    FILE* f = fopen(filename, "wb");
    assert(f != NULL);

    array1d_write(robot.bone_ids, f);
    array1d_write(robot.bone_offset_positions, f);
    array1d_write(robot.bone_rest_positions, f);
    array1d_write(robot.bone_rest_rotations, f);
    array1d_write(robot.bone_parents, f);
    
    array1d_write(robot.all_rest_positions, f);
    array1d_write(robot.all_rest_normals, f);
    array1d_write(robot.texcoords, f);
    array1d_write(robot.indices, f);
    
    array2d_write(robot.bone_weights, f);
    array2d_write(robot.bone_weights_ids, f);
    
    fclose(f);
}

void robot_forward_kinematics(
    const Robot &robot,
    slice1d<vec3> bone_anim_positions, 
    slice1d<quat> bone_anim_rotations,
    float a,
    float b,
    float c)
{
    bone_anim_positions.zero();
    bone_anim_rotations.zero();
    
    for(int i = 0; i < robot.nbones(); i++)
    {
        // Assumes bones are always sorted from root onwards
        int parent_id = robot.bone_parents(i);

        if(parent_id == -1)
        {
            bone_anim_rotations(i) = quat(1.f, 0.f, 0.f, 0.f);
            bone_anim_positions(i) = robot.bone_offset_positions(i);
        }
        else
        {
            bone_anim_rotations(i) = 
                bone_anim_rotations(parent_id) * quat(a, vec3(1.f, 0.f, 0.f)) * quat(b, vec3(0.f, 1.f, 0.f));
            bone_anim_positions(i) = 
                bone_anim_positions(parent_id) + 
                bone_anim_rotations(parent_id) * robot.bone_offset_positions(i);
        }
    }
}

Mesh make_robot_rest_mesh(const Robot &robot)
{
    Mesh mesh;

    mesh.vertex_count = robot.all_rest_positions.size;
    mesh.triangle_count = robot.indices.size / 3;
    mesh.vertices = (float *)malloc(robot.all_rest_positions.size * 3 * sizeof(float));
    mesh.normals = (float *)malloc(robot.all_rest_normals.size * 3 * sizeof(float));
    mesh.texcoords = (float *)malloc(robot.texcoords.size * 2 * sizeof(float));
    mesh.indices = (unsigned short *)malloc(robot.indices.size * sizeof(unsigned short));

    memcpy(mesh.vertices, robot.all_rest_positions.data, robot.all_rest_positions.size * 3 * sizeof(float));
    memcpy(mesh.normals, robot.all_rest_normals.data, robot.all_rest_normals.size * 3 * sizeof(float));
    memcpy(mesh.texcoords, robot.texcoords.data, robot.texcoords.size * 2 * sizeof(float));
    memcpy(mesh.indices, robot.indices.data, robot.indices.size * sizeof(unsigned short));

    return mesh;
}

void deform_robot_anim_mesh(
    const Robot &robot,
    const slice1d<vec3> bone_anim_positions,
    const slice1d<quat> bone_anim_rotations,
    Mesh &mesh)
{
    linear_blend_positions(
        slice1d<vec3>(mesh.vertex_count, (vec3 *)mesh.vertices),
        bone_anim_positions,
        bone_anim_rotations,
        robot.all_rest_positions,
        robot.bone_rest_positions,
        robot.bone_rest_rotations,
        robot.bone_weights_ids,
        robot.bone_weights);
    linear_blend_normals(
        slice1d<vec3>(mesh.vertex_count, (vec3 *)mesh.normals),
        bone_anim_rotations,
        robot.all_rest_normals,
        robot.bone_rest_rotations,
        robot.bone_weights_ids,
        robot.bone_weights);
}

Robot create_robot(float r, float h, int r_size, int h_size, int bones_size, float damping_rate = 0.8)
{
    Robot robot;

    robot.resize_bone(bones_size);

    for(int i = 0; i < bones_size; i++)
    {
        robot.bone_ids(i) = i;
        robot.bone_offset_positions(i) = vec3(0.f, h / bones_size, 0.f);
        robot.bone_rest_positions(i) = vec3(0.f, h * i / bones_size, 0.f);
        robot.bone_rest_rotations(i) = quat(1.f, 0.f, 0.f, 0.f);
        robot.bone_parents(i) = i - 1;
    }
    robot.bone_offset_positions(0) = vec3(0.f);

    int triangles_size = 2 * r_size + 2 * h_size * r_size;
    int vertices_size = 2 + (h_size + 1) * r_size;
    robot.resize_mesh(triangles_size, vertices_size);

    robot.all_rest_positions(0) = vec3(0.f, 0.f, 0.f);
    robot.all_rest_normals(0) = vec3(0.f, -1.f, 0.f);
    robot.texcoords(0) = vec2(0.f, 0.f);
    robot.bone_weights(0, 0) = 1.f;
    robot.bone_weights(0, 1) = 0.f;
    robot.bone_weights(0, 2) = 0.f;
    robot.bone_weights(0, 3) = 0.f;
    robot.bone_weights_ids(0, 0) = 0;
    robot.bone_weights_ids(0, 1) = 1;
    robot.bone_weights_ids(0, 2) = 2;
    robot.bone_weights_ids(0, 3) = 3;

    robot.all_rest_positions(vertices_size - 1) = vec3(0.f, h, 0.f);
    robot.all_rest_normals(vertices_size - 1) = vec3(0.f, 1.f, 0.f);
    robot.texcoords(vertices_size - 1) = vec2(0.f, 0.f);
    robot.bone_weights(vertices_size - 1, 0) = 1.f;
    robot.bone_weights(vertices_size - 1, 1) = 0.f;
    robot.bone_weights(vertices_size - 1, 2) = 0.f;
    robot.bone_weights(vertices_size - 1, 3) = 0.f;
    robot.bone_weights_ids(vertices_size - 1, 0) = bones_size - 1;
    robot.bone_weights_ids(vertices_size - 1, 1) = bones_size - 2;
    robot.bone_weights_ids(vertices_size - 1, 2) = bones_size - 3;
    robot.bone_weights_ids(vertices_size - 1, 3) = bones_size - 4;

    for(int i = 0; i < (h_size + 1); i++)
    {
        for(int j = 0; j < r_size; j++)
        {
            int id = 1 + i * r_size + j;
            float phi = j * 2.f * PI / r_size;
            float k = float(i) / h_size;
            robot.all_rest_positions(id) = vec3(r * sin(phi), h * k, r * cos(phi));
            robot.all_rest_normals(id) = normalize(vec3(sin(phi), 0.f, cos(phi)));
            robot.texcoords(id) = vec2(phi / 2.f / PI, k);

            std::map<int, float> temp;
            
            for(int bone_id = 0; bone_id < bones_size; bone_id++)
            {
                temp[bone_id] = length(
                                robot.all_rest_positions(id) - 
                                (robot.bone_rest_positions(bone_id) + robot.bone_offset_positions(bone_id)));

                // std::cout << temp[bone_id] << ',';
            }

            // std::cout << std::endl;

            std::vector<std::pair<int, float>> distances(temp.begin(), temp.end());

            sort(distances.begin(), distances.end(), compare_map_int_float);

            // std::cout << distances[0].second << ',' << distances[1].second << ',' 
            //           << distances[2].second << ',' << distances[3].second << std::endl;
            // std::cout << distances[0].first << ',' << distances[1].first << ',' 
            //           << distances[2].first << ',' << distances[3].first << std::endl;
            
            vec4 bone_weight = to_one(vec4(
                                        1.f / distances[0].second, 
                                        1.f / distances[1].second, 
                                        1.f / distances[2].second, 
                                        1.f / distances[3].second));

            // vec4 bone_weight = normalize(vec4(
            //                             1.f, 
            //                             0.f, 
            //                             0.f, 
            //                             0.f));

            robot.bone_weights(id, 0) = bone_weight.x;
            robot.bone_weights(id, 1) = bone_weight.y;
            robot.bone_weights(id, 2) = bone_weight.z;
            robot.bone_weights(id, 3) = bone_weight.w;
            robot.bone_weights_ids(id, 0) = distances[0].first;
            robot.bone_weights_ids(id, 1) = distances[1].first;
            robot.bone_weights_ids(id, 2) = distances[2].first;
            robot.bone_weights_ids(id, 3) = distances[3].first;
        }
    }

    for(int i = 0; i < r_size; i++)
    {
        robot.indices(3 * i + 0) = 0;
        robot.indices(3 * i + 1) = i + 2;
        robot.indices(3 * i + 2) = i + 1;
    }

    // printf("%d, %d\n", 0, 3 * (r_size - 1) + 2);

    for(int i = 0; i < r_size; i++)
    {
        int id = r_size + 2 * h_size * r_size + i;
        robot.indices(3 * id + 0) = 1 + (h_size + 1) * r_size;
        robot.indices(3 * id + 1) = h_size * r_size + circulate_int(i + 1, 1, r_size);
        robot.indices(3 * id + 2) = h_size * r_size + circulate_int(i + 2, 1, r_size);
    }

    // printf("%d, %d\n", 3 * (r_size + 2 * h_size * r_size), 3 * (r_size + 2 * h_size * r_size + r_size - 1) + 2);

    for(int i = 0; i < h_size; i++)
    {
        for(int j = 0; j < r_size; j++)
        {
            int id = r_size + 2 * i * r_size + 2 * j;
            
            int v0 = (i + 1) * r_size + circulate_int(j + 1, 1, r_size);
            int v1 = i * r_size + circulate_int(j + 1, 1, r_size);
            int v2 = i * r_size + circulate_int(j + 2, 1, r_size);
            int v3 = (i + 1) * r_size + circulate_int(j + 2, 1, r_size);

            robot.indices(3 * id + 0) = v0;
            robot.indices(3 * id + 1) = v1;
            robot.indices(3 * id + 2) = v2;
            robot.indices(3 * id + 3) = v2;
            robot.indices(3 * id + 4) = v3;
            robot.indices(3 * id + 5) = v0;
        }
    }

    // printf("%d, %d\n", 3 * r_size, 3 * (r_size + 2 * (h_size - 1) * r_size + 2 * (r_size - 1)) + 5);

    return robot;
}

#endif