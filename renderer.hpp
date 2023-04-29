#ifndef ENGINE_RENDERER
#define ENGINE_RENDERER

#include "camera.hpp"
#include "light.hpp"
#include "model.hpp"
#include "mesh.hpp"

enum G_buffer
{
    color_buffer       = 0,
    vertice_buffer     = 1,
    normal_buffer      = 2,
    depth_buffer       = 3,
};

struct Renderer
{
    std::vector<Model *> models;
    std::vector<Point_Light *> point_lights;
    std::vector<Polygon_Light *> polygon_lights;
};

struct updated_paramters
{
    vec3 *light_pos;
    vec3 *light_dir;
    vec3 *light_radiance;

    FBO *shadow_map;
};

void bind_geometry_info(
    Shader &shader, 
    int triangle_count, 
    unsigned short *indices, 
    vec3 *vertices, 
    vec3 *normals, 
    vec2 *texcoords)
{
    for(int i = 0; i < triangle_count; i++)
    {
        // std::cout << i << std::endl;
        int t0 = indices[3 * i];
        int t1 = indices[3 * i + 1];
        int t2 = indices[3 * i + 2];

        // printf("%d, %d, %d\n", t0, t1, t2);

        // print(vertices[t0]);
        // print(vertices[t1]);
        // print(vertices[t2]);

        shader.ver[i].vertex_pos[0] = vertices[t0];
        shader.ver[i].vertex_pos[1] = vertices[t1];
        shader.ver[i].vertex_pos[2] = vertices[t2];

        shader.ver[i].vertex_normal[0] = normals[t0];
        shader.ver[i].vertex_normal[1] = normals[t1];
        shader.ver[i].vertex_normal[2] = normals[t2];

        shader.ver[i].vertex_texcoords[0] = texcoords[t0];
        shader.ver[i].vertex_texcoords[1] = texcoords[t1];
        shader.ver[i].vertex_texcoords[2] = texcoords[t2];
    }

    // std::cout << "geo!" << std::endl;
}

void bind_camera_parameters(Shader &shader, Model *model, Camera &camera)
{
    shader.model_matrix = model->transform;
    shader.view_matrix = camera.get_view_matrix();
    shader.projection_matrix = camera.get_projection_matrix();

    shader.camera_pos = &camera.position;
    shader.camera_dir = &camera.direction;
    shader.G_buffer = camera.fbo;

    // print(camera.position);
    // print(camera.direction);

    // std::cout << "cam" << std::endl;
}

void bind_material_parameters(Shader &shader, Material &material)
{
    shader.albedo_map = material.tex;
    shader.kd = material.kd;
    shader.ks = material.ks;
    shader.ka = material.ka;
    shader.metallic = material.metallic;
    shader.roughness = material.roughness;

    // std::cout << "mat" << std::endl;
}

void update_material_parameters(Shader &shader, updated_paramters *par)
{
    shader.shadow_map = par->shadow_map;
    shader.light_pos = par->light_pos;
    shader.light_dir = par->light_dir;
    shader.light_radiance = par->light_radiance;

    // std::cout << "update" << std::endl;
}

void rasterizer(Shader &shader, FBO *fbo)
{
    FBO depth_map(fbo->rows, fbo->cols, vec3(100.f));

    mat4 modelview_matrix = shader.view_matrix * shader.model_matrix;
    mat4 projection_matrix = shader.projection_matrix;

    // print(modelview_matrix);
    // print(projection_matrix);

    for(int i = 0; i < shader.triangle_count; i++)
    {
        // std::cout << i << std::endl;

        vec4 v0 = modelview_matrix * pos_to_vec4(shader.ver[i].vertex_pos[0]);
        vec4 v1 = modelview_matrix * pos_to_vec4(shader.ver[i].vertex_pos[1]);
        vec4 v2 = modelview_matrix * pos_to_vec4(shader.ver[i].vertex_pos[2]);

        // print(v0);
        // print(v1);
        // print(v2);

        // print(projection_matrix * v0);
        // print(projection_matrix * v1);
        // print(projection_matrix * v2);

        vec4 p0 = normalize(projection_matrix * v0);
        vec4 p1 = normalize(projection_matrix * v1);
        vec4 p2 = normalize(projection_matrix * v2);

        // print(p0);
        // print(p1);
        // print(p2);

        float max_x = clampf(0.5f * fbo->cols * (1.f + maxf_3(p0.x, p1.x, p2.x)), 0, fbo->cols - 1);
        float min_x = clampf(0.5f * fbo->cols * (1.f + minf_3(p0.x, p1.x, p2.x)), 0, fbo->cols - 1);
        float max_y = clampf(0.5f * fbo->rows * (1.f + maxf_3(p0.y, p1.y, p2.y)), 0, fbo->rows - 1);
        float min_y = clampf(0.5f * fbo->rows * (1.f + minf_3(p0.y, p1.y, p2.y)), 0, fbo->rows - 1);

        // printf("\n%f, %f, %f, %f\n", max_x, min_x, max_y, min_y);

        for(int x = min_x; x < max_x; x++)
        {
            for(int y = min_y; y < max_y; y++)
            {
                int fbo_y = fbo->rows - 1 - y;
                int fbo_x = x;

                vec2 p = vec2(x * 2.f / fbo->cols - 1.f, y * 2.f / fbo->rows - 1.f);

                // print(p);

                float alpha, beta, gamma;
                bool flag = inside_2dtriangle(vec4_to_vec2(p0), vec4_to_vec2(p1), vec4_to_vec2(p2), p, alpha, beta, gamma);
                // std::cout << flag << std::endl;
                float depth = -(alpha * v0.z + beta * v1.z + gamma * v2.z);

                if(flag == true && depth < depth_map(fbo_y, fbo_x).z)
                {
                    fragment_payload frag;

                    vec3 color;

                    depth_map(fbo_y, fbo_x).z = depth;
                    frag.depth = depth;
                    
                    update_fragment_payload(shader.ver[i], frag, alpha, beta, gamma);
                    color = 255.f * shader.pbr_shader(frag);

                    fbo->colors(fbo_y, fbo_x) = clampv(color);   

                    // print(fbo->colors(fbo_y, fbo_x));   
                }
            }
        }
    }
}

void draw(Model *model, Camera &camera, FBO *fbo, updated_paramters *par)
{
    for(int i = 0; i < model->mesh_count; i++)
    {
        Shader shader(model->meshes[i].triangle_count);

        // std::cout << shader.triangle_count << std::endl;

        bind_geometry_info(
            shader, 
            model->meshes[i].triangle_count, 
            model->meshes[i].indices,
            (vec3 *)model->meshes[i].vertices,
            (vec3 *)model->meshes[i].normals,
            (vec2 *)model->meshes[i].texcoords);
        bind_camera_parameters(
            shader, 
            model, 
            camera);
        bind_material_parameters(
            shader, 
            model->materials[i]);
        update_material_parameters(
            shader, 
            par);
        rasterizer(
            shader, 
            fbo);
    }
}

#endif