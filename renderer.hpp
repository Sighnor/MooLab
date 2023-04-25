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
    std::vector<Model> models;
    std::vector<Point_Light> point_lights;
    std::vector<Polygon_Light> polygon_lights;
};

struct updated_paramters
{
    vec3 light_pos;
    vec3 light_dir;
    vec3 light_radiance;

    FBO *shadow_map;
};

void bind_geometry_info(
    Shader &shader, 
    int triangle_count, 
    int *indices, 
    vec3 *vertices, 
    vec3 *normals, 
    vec2 *texcoords)
{
    for(int i = 0; i < triangle_count; i++)
    {
        int t0 = indices[3 * i];
        int t1 = indices[3 * i + 1];
        int t2 = indices[3 * i + 2];

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
}

void bind_camera_parameters(Shader &shader, Model &model, Camera &camera)
{
    shader.model_matrix = model.transform;
    shader.view_matrix = camera.get_view_matrix();
    shader.projection_matrix = camera.get_projection_matrix();

    shader.camera_pos = camera.position;
    shader.camera_dir = camera.direction;
    shader.G_buffer = camera.fbo;
}

void bind_material_parameters(Shader &shader, Material &material)
{
    shader.albedo_map = material.tex;
    shader.kd = material.kd;
    shader.ks = material.ks;
    shader.ka = material.ka;
    shader.metallic = material.metallic;
    shader.roughness = material.roughness;
}

void update_material_parameters(Shader &shader, updated_paramters *par)
{
    shader.shadow_map = par->shadow_map;
    shader.light_pos = par->light_pos;
    shader.light_dir = par->light_dir;
    shader.light_radiance = par->light_radiance;
}

void rasterizer(Shader shader, FBO *fbo)
{
    FBO depth_map;

    mat4 mvp = shader.projection_matrix * shader.view_matrix * shader.model_matrix;

    for(int i = 0; i < shader.triangle_count; i++)
    {
        vec4 v0 = normalize(mvp * pos_to_vec4(shader.ver[i].vertex_pos[0]));
        vec4 v1 = normalize(mvp * pos_to_vec4(shader.ver[i].vertex_pos[1]));
        vec4 v2 = normalize(mvp * pos_to_vec4(shader.ver[i].vertex_pos[2]));

        float max_x = maxf_3(v0.x, v1.x, v2.x);
        float min_x = minf_3(v0.x, v1.x, v2.x);
        float max_y = maxf_3(v0.y, v1.y, v2.y);
        float min_y = minf_3(v0.y, v1.y, v2.y);

        for(int x = min_x; x < max_x; x++)
        {
            for(int y = min_y; y < max_y; y++)
            {
                float alpha, beta, gamma;
                float depth = inside_triangle(vec4_to_vec2(v0), vec4_to_vec2(v0), vec4_to_vec2(v0), vec2(x, y), alpha, beta, gamma);

                if(depth > depth_map(x, y).z)
                {
                    fragment_payload frag;

                    vec3 color;

                    depth_map(x, y).z = depth;
                    frag.depth = depth;
                    
                    update_fragment_payload(shader.ver[i], frag, alpha, beta, gamma);
                    color = pbr_shader(shader, frag);

                    fbo->colors(x, y) = clamp(fbo->colors(x, y) + color);      
                }
            }
        }
    }
}

void draw(Model model, Camera camera, FBO *fbo, updated_paramters *par)
{
    for(int i = 0; i < model.mesh_count; i++)
    {
        Shader shader(model.meshes[i].triangle_count);

        bind_geometry_info(
            shader, 
            model.meshes[i].triangle_count, 
            model.meshes[i].indices,
            (vec3 *)model.meshes[i].vertices,
            (vec3 *)model.meshes[i].vertices,
            (vec2 *)model.meshes[i].texcoords);
        bind_camera_parameters(
            shader, 
            model, 
            camera);
        bind_material_parameters(
            shader, 
            model.materials[i]);
        update_material_parameters(
            shader, 
            par);
        rasterizer(
            shader, 
            fbo);
    }
}

#endif