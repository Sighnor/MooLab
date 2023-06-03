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
    std::vector<Model*> models;
    std::vector<Point_Light*> point_lights;
    std::vector<Polygon_Light*> polygon_lights;
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

void bind_camera_parameters(Shader &shader, Model *model, Camera &camera)
{
    shader.model_matrix = model->transform;
    shader.view_matrix = camera.get_view_matrix();
    shader.projection_matrix = camera.get_projection_matrix();
    shader.screen_to_world_matrix = camera.get_inverse_view_matrix() * camera.get_inverse_projection_matrix();

    shader.camera_pos = camera.position;
    shader.camera_dir = camera.direction;
    shader.z_near = camera.z_near;
    shader.G_buffer = camera.fbo;
}

void bind_material_parameters(Shader &shader, Material &material)
{
    shader.albedo_map = &material.tex;
    shader.BRDFLut = &material.BRDFLut;
    shader.EavgLut = &material.EavgLut;
    shader.kd = material.kd;
    shader.ks = material.ks;
    shader.ka = material.ka;
    shader.metallic = material.metallic;
    shader.roughness = material.roughness;
}

void update_material_parameters(Shader &shader, updated_paramters *par)
{
    shader.shadow_map = par->shadow_map;
    shader.light_pos = *par->light_pos;
    shader.light_dir = *par->light_dir;
    shader.light_radiance = *par->light_radiance;
}

void rasterizer(Shader &shader, FBO *fbo, bool ifdraw)
{
    mat4 modelview_matrix = shader.view_matrix * shader.model_matrix;
    mat4 projection_matrix = shader.projection_matrix;

    #pragma omp parallel for 
    for(int i = 0; i < shader.triangle_count; i++)
    {
        vec4 v0 = modelview_matrix * pos_to_vec4(shader.ver[i].vertex_pos[0]);
        vec4 v1 = modelview_matrix * pos_to_vec4(shader.ver[i].vertex_pos[1]);
        vec4 v2 = modelview_matrix * pos_to_vec4(shader.ver[i].vertex_pos[2]);

        if(- v0.z < shader.z_near && - v1.z < shader.z_near && - v2.z < shader.z_near)
        {
            continue;
        }

        vec4 p0 = standardize(projection_matrix * v0);
        vec4 p1 = standardize(projection_matrix * v1);
        vec4 p2 = standardize(projection_matrix * v2);

        float max_x = clampf(0.5f * fbo->cols * (1.f + maxf_3(p0.x, p1.x, p2.x)), 0, fbo->cols - 1);
        float min_x = clampf(0.5f * fbo->cols * (1.f + minf_3(p0.x, p1.x, p2.x)), 0, fbo->cols - 1);
        float max_y = clampf(0.5f * fbo->rows * (1.f + maxf_3(p0.y, p1.y, p2.y)), 0, fbo->rows - 1);
        float min_y = clampf(0.5f * fbo->rows * (1.f + minf_3(p0.y, p1.y, p2.y)), 0, fbo->rows - 1);

        for(int x = min_x; x < max_x; x++)
        {
            for(int y = min_y; y < max_y; y++)
            {
                int fbo_y = fbo->rows - 1 - (y + 0.5f);
                int fbo_x = x + 0.5f;

                vec2 p = vec2((x + 0.5f) * 2.f / fbo->cols - 1.f, (y + 0.5f) * 2.f / fbo->rows - 1.f);

                float alpha, beta, gamma;
                bool flag = inside_2dtriangle(vec4_to_vec2(p0), vec4_to_vec2(p1), vec4_to_vec2(p2), p, alpha, beta, gamma);

                float depth = -(alpha * v0.z + beta * v1.z + gamma * v2.z);

                if(flag == true && depth <= fbo->getdepth(fbo_y, fbo_x) && depth > shader.z_near)
                {
                    fbo->getdepth(fbo_y, fbo_x) = depth;

                    if(ifdraw == true)
                    {
                        switch (shader.shader_type)
                        {
                            case Light:
                            {
                                fbo->getcolor(fbo_y, fbo_x) = clampv(255.f * shader.light_radiance);
                                break;
                            }
                            case PBR:
                            {   
                                fragment_payload frag = shader.get_fragment_payload(i, alpha, beta, gamma, depth);

                                vec3 color;
                                color = 255.f * shader.pbr_shader(frag);
                                fbo->getcolor(fbo_y, fbo_x) = clampv(color);
                                // fbo->getcolor(fbo_y, fbo_x) = clampv(255.f * 0.5f);
                                // fbo->getcolor(fbo_y, fbo_x) = 255.f * fbo->getdepth(fbo_y, fbo_x);
                                break;
                            }
                            default:
                                break;
                        }
                    }
                }
            }
        }
    }
}

void draw(Model *model, Camera &camera, FBO *fbo, updated_paramters *par, bool ifdraw)
{
    for(int i = 0; i < model->mesh_count; i++)
    {
        bind_geometry_info(
            model->shaders[i], 
            model->meshes[i].triangle_count, 
            model->meshes[i].indices,
            (vec3 *)model->meshes[i].vertices,
            (vec3 *)model->meshes[i].normals,
            (vec2 *)model->meshes[i].texcoords);
        bind_camera_parameters(
            model->shaders[i], 
            model, 
            camera);
        bind_material_parameters(
            model->shaders[i], 
            model->materials[i]);
        update_material_parameters(
            model->shaders[i], 
            par);
        rasterizer(
            model->shaders[i], 
            fbo,
            ifdraw);
    }
}

void draw_point(vec3 point, Camera &camera, FBO *fbo, vec3 color)
{
    mat4 view_matrix = camera.get_view_matrix();
    mat4 projection_matrix = camera.get_projection_matrix();
    vec4 p = standardize(projection_matrix * view_matrix * pos_to_vec4(point));

    int x = 0.5f * fbo->cols * (1.f + p.x);
    int y = 0.5f * fbo->rows * (1.f + p.y);

    int fbo_x = x;
    int fbo_y = fbo->rows - 1 - y;

    if(fbo_y >= 0 && fbo_y <= fbo->rows - 1 && fbo_x >= 0 && fbo_x <= fbo->cols - 1)
    {
        for(int i = std::max(fbo_x - 2, 0); i < std::min(fbo_x + 3, fbo->cols - 1); i++)
        {
            for(int j = std::max(fbo_y - 2, 0); j < std::min(fbo_y + 3, fbo->rows - 1); j++)
            {
                fbo->getcolor(j, i) = clampv(color);
            }
        }
    }
}

#endif