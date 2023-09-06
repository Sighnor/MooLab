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

struct MooRenderer
{
    std::vector<MooModel*> models;
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
    MooShader &shader, 
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

void bind_camera_parameters(MooShader &shader, MooModel *model, MooCamera &camera)
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

void bind_material_parameters(MooShader &shader, MooMaterial &material)
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

void update_material_parameters(MooShader &shader, updated_paramters *par)
{
    shader.shadow_map = par->shadow_map;
    shader.light_pos = *par->light_pos;
    shader.light_dir = *par->light_dir;
    shader.light_radiance = *par->light_radiance;
}

bool inside_2dtriangle(const vec2 &p0, const vec2 &p1, const vec2 &p2, const vec2 &p, float &alpha, float &beta, float &gamma)
{
    alpha = ((p1.y - p2.y) * p.x + (p2.x - p1.x) * p.y + p1.x * p2.y - p2.x * p1.y) / ((p1.y - p2.y) * p0.x + (p2.x - p1.x) * p0.y + p1.x * p2.y - p2.x * p1.y);
    beta = ((p2.y - p0.y) * p.x + (p0.x - p2.x) * p.y + p2.x * p0.y - p0.x * p2.y) / ((p2.y - p0.y) * p1.x + (p0.x - p2.x) * p1.y + p2.x * p0.y - p0.x * p2.y);
    gamma = ((p0.y - p1.y) * p.x + (p1.x - p0.x) * p.y + p0.x * p1.y - p1.x * p0.y) / ((p0.y - p1.y) * p2.x + (p1.x - p0.x) * p2.y + p0.x * p1.y - p1.x * p0.y);

    if(alpha >= 0.f && alpha <= 1.f && beta >= 0.f && beta <= 1.f && gamma >= 0.f && gamma <= 1.f)
    {
        return true;
    }

    return false;
}

bool inside_3dtriangle(const vec3 &v0, const vec3 &v1, const vec3 &v2, const vec3 &v, float &alpha, float &beta, float &gamma)
{
    //v0, v1, v2逆时针
    mat3 tri_matrix = inv_mat(get_coordinate_matrix(cross(v1 - v0, v2 - v0), vec3(0.f, 1.f, 0.f)));

    vec2 p = vec3_to_vec2(tri_matrix * v);
    vec2 p0 = vec3_to_vec2(tri_matrix * v0);
    vec2 p1 = vec3_to_vec2(tri_matrix * v1);
    vec2 p2 = vec3_to_vec2(tri_matrix * v2);

    return inside_2dtriangle(p0, p1, p2, p, alpha, beta, gamma);
}

void rasterizer(MooShader &shader, FBO *fbo, bool ifdraw)
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
                            case LIGHT:
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
                            case TEXTURE:
                            {
                                fragment_payload frag = shader.get_fragment_payload(i, alpha, beta, gamma, depth);

                                vec3 color;
                                color = 255.f * shader.texture_shader(frag);
                                fbo->getcolor(fbo_y, fbo_x) = clampv(color);
                            }
                            case NORMALS:
                            {
                                fragment_payload frag = shader.get_fragment_payload(i, alpha, beta, gamma, depth);

                                vec3 color;
                                color = 255.f * shader.normal_shader(frag);
                                fbo->getcolor(fbo_y, fbo_x) = clampv(color);
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

void draw(MooModel *model, MooCamera &camera, FBO *fbo, updated_paramters *par, bool ifdraw)
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

void draw_point(vec3 point, MooCamera &camera, FBO *fbo, vec3 color)
{
    mat4 view_matrix = camera.get_view_matrix();
    mat4 projection_matrix = camera.get_projection_matrix();
    vec4 v = view_matrix * pos_to_vec4(point);
    if(- v.z < camera.z_near)
    {
        return;
    }
    vec4 p = standardize(projection_matrix * view_matrix * pos_to_vec4(point));

    int x = 0.5f * fbo->cols * (1.f + p.x);
    int y = fbo->rows - 1 - 0.5f * fbo->rows * (1.f + p.y);

    if(y >= 0 && y <= fbo->rows - 1 && x >= 0 && x <= fbo->cols - 1)
    {
        for(int i = std::max(x - 2, 0); i < std::min(x + 3, fbo->cols - 1); i++)
        {
            for(int j = std::max(y - 2, 0); j < std::min(y + 3, fbo->rows - 1); j++)
            {
                fbo->getcolor(j, i) = clampv(color);
            }
        }
    }
}

// Bresenham's line drawing algorithm
// Code taken from a stack overflow answer: https://stackoverflow.com/a/16405254
void draw_line(vec3 begin, vec3 end, MooCamera &camera, FBO *fbo, vec3 color)
{
    mat4 view_matrix = camera.get_view_matrix();
    mat4 projection_matrix = camera.get_projection_matrix();
    vec4 v0 = view_matrix * pos_to_vec4(begin);
    vec4 v1 = view_matrix * pos_to_vec4(end);
    if(- v0.z < camera.z_near || - v1.z < camera.z_near)
    {
        return;
    }
    vec4 p0 = standardize(projection_matrix * view_matrix * pos_to_vec4(begin));
    vec4 p1 = standardize(projection_matrix * view_matrix * pos_to_vec4(end));

    float x0 = 0.5f * fbo->cols * (1.f + p0.x);
    float y0 = fbo->rows - 1 - 0.5f * fbo->rows * (1.f + p0.y);
    float x1 = 0.5f * fbo->cols * (1.f + p1.x);
    float y1 = fbo->rows - 1 - 0.5f * fbo->rows * (1.f + p1.y);

    int x, y, dx, dy, dx0, dy0, px, py, xe, ye, i;

    dx = x1 - x0;
    dy = y1 - y0;
    dx0 = fabs(dx);
    dy0 = fabs(dy);
    px = 2 * dy0 - dx0;
    py = 2 * dx0 - dy0;

    if(dy0 <= dx0)
    {
        if(dx >= 0)
        {
            x = x0;
            y = y0;
            xe = x1;
        }
        else
        {
            x = x1;
            y = y1;
            xe = x0;
        }
        fbo->getcolor(y, x) = clampv(color);
        for(i = 0; x < xe; i++)
        {
            x = x + 1;
            if(px < 0)
            {
                px = px + 2 * dy0;
            }
            else
            {
                if((dx < 0 && dy < 0) || (dx > 0 && dy > 0))
                {
                    y = y + 1;
                }
                else
                {
                    y = y - 1;
                }
                px = px + 2 * (dy0 - dx0);
            }
            fbo->getcolor(y, x) = clampv(color);
        }
    }
    else
    {
        if(dy >= 0)
        {
            x = x0;
            y = y0;
            ye = y1;
        }
        else
        {
            x = x1;
            y = y1;
            ye = y0;
        }
        fbo->getcolor(y, x) = clampv(color);
        for(i = 0; y < ye; i++)
        {
            y = y + 1;
            if(py <= 0)
            {
                py = py + 2 * dx0;
            }
            else
            {
                if((dx < 0 && dy < 0) || (dx > 0 && dy > 0))
                {
                    x = x + 1;
                }
                else
                {
                    x = x - 1;
                }
                py = py + 2 * (dx0 - dy0);
            }
            fbo->getcolor(y, x) = clampv(color);
        }
    }
}

#endif