#ifndef MOOLAB_RENDERER
#define MOOLAB_RENDERER

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

//三角形建立阶段的结果，光栅阶段只读
struct FRasterTriangle
{
    vec2 P0, P1, P2;                    //NDC坐标xy
    float Z0, Z1, Z2;                   //观察空间z，用于插值深度
    float MinX, MaxX, MinY, MaxY;       //像素包围盒，边界与比较方式与原实现一致
    bool bVisible;
};

//每个行带的像素行数；行带之间写入不同的FBO行，可无锁并行
constexpr int RasterBandRows = 8;

//两阶段光栅：
//  1. 三角形建立（按三角形并行）：变换、近平面剔除、包围盒；随后串行分箱到各行带
//  2. 光栅+着色（按屏幕行带并行）：每个行带按原三角形顺序遍历本箱，只写本行带的像素
//原实现按三角形并行并直接读写深度/颜色，重叠三角形间存在数据竞争；
//现实现中每个像素仍按三角形序号依次做深度测试，结果与串行执行逐位一致且与线程数无关。
//使用函数内静态缓冲，仅支持单一渲染线程调用（不可重入）。shader以const传入，由类型保证并行区内只读。
void rasterizer(const MooShader &shader, FBO *fbo, bool ifdraw)
{
    mat4 modelview_matrix = shader.view_matrix * shader.model_matrix;
    mat4 projection_matrix = shader.projection_matrix;

    const int triangle_count = shader.triangle_count;
    const int band_count = (fbo->rows + RasterBandRows - 1) / RasterBandRows;

    static std::vector<FRasterTriangle> triangles;
    static std::vector<std::vector<int>> band_bins;
    triangles.resize(std::max(triangle_count, 0));
    band_bins.resize(band_count);

    #pragma omp parallel
    {
        #pragma omp for
        for(int i = 0; i < triangle_count; i++)
        {
            FRasterTriangle &tri = triangles[i];

            vec4 v0 = modelview_matrix * pos_to_vec4(shader.ver[i].vertex_pos[0]);
            vec4 v1 = modelview_matrix * pos_to_vec4(shader.ver[i].vertex_pos[1]);
            vec4 v2 = modelview_matrix * pos_to_vec4(shader.ver[i].vertex_pos[2]);

            if(- v0.z < shader.z_near && - v1.z < shader.z_near && - v2.z < shader.z_near)
            {
                tri.bVisible = false;
                continue;
            }

            vec4 p0 = standardize(projection_matrix * v0);
            vec4 p1 = standardize(projection_matrix * v1);
            vec4 p2 = standardize(projection_matrix * v2);

            tri.P0 = vec4_to_vec2(p0);
            tri.P1 = vec4_to_vec2(p1);
            tri.P2 = vec4_to_vec2(p2);
            tri.Z0 = v0.z;
            tri.Z1 = v1.z;
            tri.Z2 = v2.z;
            tri.MaxX = clampf(0.5f * fbo->cols * (1.f + maxf_3(p0.x, p1.x, p2.x)), 0, fbo->cols - 1);
            tri.MinX = clampf(0.5f * fbo->cols * (1.f + minf_3(p0.x, p1.x, p2.x)), 0, fbo->cols - 1);
            tri.MaxY = clampf(0.5f * fbo->rows * (1.f + maxf_3(p0.y, p1.y, p2.y)), 0, fbo->rows - 1);
            tri.MinY = clampf(0.5f * fbo->rows * (1.f + minf_3(p0.y, p1.y, p2.y)), 0, fbo->rows - 1);
            //w接近0时包围盒可能为NaN，转int属UB且会导致分箱越界；原实现对此类三角形不会画出有效像素，直接剔除
            tri.bVisible = std::isfinite(tri.MinX) && std::isfinite(tri.MaxX) && std::isfinite(tri.MinY) && std::isfinite(tri.MaxY);
        }
        //omp for结束处的隐式屏障保证建立阶段全部完成

        //串行分箱，箱内三角形序号递增（保证深度测试顺序与串行一致）；行范围可偏保守，精确边界在光栅阶段判定
        #pragma omp single
        {
            for(std::vector<int> &bin : band_bins)
            {
                bin.clear();
            }
            for(int i = 0; i < triangle_count; i++)
            {
                const FRasterTriangle &tri = triangles[i];
                if(!tri.bVisible)
                {
                    continue;
                }
                const int y_first = static_cast<int>(tri.MinY);
                if(!(y_first < tri.MaxY))
                {
                    continue;
                }
                const int band_first = y_first / RasterBandRows;
                const int band_last = std::min(static_cast<int>(tri.MaxY) / RasterBandRows, band_count - 1);
                MOO_CHECK(band_first >= 0 && band_first <= band_last);
                for(int band = band_first; band <= band_last; band++)
                {
                    band_bins[band].push_back(i);
                }
            }
        }
        //omp single结束处的隐式屏障保证分箱完成

        #pragma omp for schedule(dynamic, 1)
        for(int band = 0; band < band_count; band++)
        {
            const int band_begin = band * RasterBandRows;
            const int band_end = std::min(band_begin + RasterBandRows, fbo->rows);

            for(const int i : band_bins[band])
            {
                const FRasterTriangle &tri = triangles[i];

                //y的取值范围与原实现的 for(int y = min_y; y < max_y; y++) 一致，再与本行带求交
                const int y_begin = std::max(static_cast<int>(tri.MinY), band_begin);
                if(y_begin >= band_end || !(y_begin < tri.MaxY))
                {
                    continue;
                }

                for(int y = y_begin; y < tri.MaxY && y < band_end; y++)
                {
                    for(int x = static_cast<int>(tri.MinX); x < tri.MaxX; x++)
                    {
                        int fbo_y = fbo->rows - 1 - (y + 0.5f);
                        int fbo_x = x + 0.5f;

                        vec2 p = vec2((x + 0.5f) * 2.f / fbo->cols - 1.f, (y + 0.5f) * 2.f / fbo->rows - 1.f);

                        float alpha, beta, gamma;
                        bool flag = inside_2dtriangle(tri.P0, tri.P1, tri.P2, p, alpha, beta, gamma);

                        float depth = -(alpha * tri.Z0 + beta * tri.Z1 + gamma * tri.Z2);

                        if(flag == true && depth <= fbo->getdepth(fbo_y, fbo_x) && depth > shader.z_near)
                        {
                            fbo->getdepth(fbo_y, fbo_x) = depth;

                            if(ifdraw == true)
                            {
                                switch (shader.shader_type)
                                {
                                    case LIGHT:
                                    {
                                        vec3 color;
                                        color = 255.f * shader.light_radiance;
                                        fbo->getcolor(fbo_y, fbo_x) = clampv(color);
                                        break;
                                    }
                                    case PBR:
                                    {
                                        fragment_payload frag = shader.get_fragment_payload(i, alpha, beta, gamma, depth);

                                        vec3 color;
                                        color = 255.f * shader.pbr_shader(frag);
                                        fbo->getcolor(fbo_y, fbo_x) = clampv(color);
                                        break;
                                    }
                                    case TEXTURE:
                                    {
                                        fragment_payload frag = shader.get_fragment_payload(i, alpha, beta, gamma, depth);

                                        vec3 color;
                                        color = 255.f * shader.texture_shader(frag);
                                        fbo->getcolor(fbo_y, fbo_x) = clampv(color);
                                        break;
                                    }
                                    case NORMALS:
                                    {
                                        fragment_payload frag = shader.get_fragment_payload(i, alpha, beta, gamma, depth);

                                        vec3 color;
                                        color = 255.f * shader.normal_shader(frag);
                                        fbo->getcolor(fbo_y, fbo_x) = clampv(color);
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

//draw_point的投影部分：点在近平面之前或落在画面外时返回false；(x, y)为像素坐标，y行0为顶部。
//GPU路径用它算出与CPU完全相同的像素位置后再提交绘制
inline bool ProjectPointToPixel(vec3 Point, MooCamera &Camera, int Rows, int Cols, int &OutX, int &OutY)
{
    const mat4 ViewMatrix = Camera.get_view_matrix();
    const mat4 ProjectionMatrix = Camera.get_projection_matrix();
    const vec4 ViewPos = ViewMatrix * pos_to_vec4(Point);
    if(- ViewPos.z < Camera.z_near)
    {
        return false;
    }
    const vec4 NdcPos = standardize(ProjectionMatrix * ViewMatrix * pos_to_vec4(Point));

    // 整个表达式向零截断，保持旧draw_point取整语义。
    OutX = static_cast<int>(0.5f * Cols * (1.f + NdcPos.x));
    OutY = static_cast<int>(Rows - 1 - 0.5f * Rows * (1.f + NdcPos.y));

    return OutY >= 0 && OutY <= Rows - 1 && OutX >= 0 && OutX <= Cols - 1;
}

void draw_point(vec3 point, MooCamera &camera, FBO *fbo, vec3 color, int size = 2)
{
    int x = 0;
    int y = 0;
    if(ProjectPointToPixel(point, camera, fbo->rows, fbo->cols, x, y))
    {
        for(int i = std::max(x - size, 0); i <= std::min(x + size, fbo->cols - 1); i++)
        {
            for(int j = std::max(y - size, 0); j <= std::min(y + size, fbo->rows - 1); j++)
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

    int x0 = clampf(0.5f * fbo->cols * (1.f + p0.x), 0, fbo->cols - 1);
    int y0 = clampf(fbo->rows - 1 - 0.5f * fbo->rows * (1.f + p0.y), 0, fbo->rows - 1);
    int x1 = clampf(0.5f * fbo->cols * (1.f + p1.x), 0, fbo->cols - 1);
    int y1 = clampf(fbo->rows - 1 - 0.5f * fbo->rows * (1.f + p1.y), 0, fbo->rows - 1);

    fbo->getcolor(y0, x0) = clampv(color);
    fbo->getcolor(y1, x1) = clampv(color);

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

void draw_curve(vec3 pos0, vec3 vel0, vec3 vel5, MooCamera &camera, FBO *fbo, vec3 color, int N, float dt = 0.0166667f)
{
    vec3 pos1 = pos0 + ((vel5 - vel0) / 2.f * 0.2f * 0.2f + vel0 * 0.2f) * N * dt;
    vec3 pos2 = pos0 + ((vel5 - vel0) / 2.f * 0.4f * 0.4f + vel0 * 0.4f) * N * dt;
    vec3 pos3 = pos0 + ((vel5 - vel0) / 2.f * 0.6f * 0.6f + vel0 * 0.6f) * N * dt;
    vec3 pos4 = pos0 + ((vel5 - vel0) / 2.f * 0.8f * 0.8f + vel0 * 0.8f) * N * dt;
    vec3 pos5 = pos0 + ((vel5 - vel0) / 2.f * 1.0f * 1.0f + vel0 * 1.0f) * N * dt;

    // mat4 view_matrix = camera.get_view_matrix();
    // mat4 projection_matrix = camera.get_projection_matrix();
    // vec4 v0 = view_matrix * pos_to_vec4(pos0);
    // vec4 v1 = view_matrix * pos_to_vec4(pos1);
    // vec4 v2 = view_matrix * pos_to_vec4(pos2);
    // vec4 v3 = view_matrix * pos_to_vec4(pos3);
    // vec4 v4 = view_matrix * pos_to_vec4(pos4);
    // vec4 v5 = view_matrix * pos_to_vec4(pos5);

    // if(- v0.z < camera.z_near || - v1.z < camera.z_near ||
    //    - v2.z < camera.z_near || - v3.z < camera.z_near ||
    //    - v4.z < camera.z_near || - v5.z < camera.z_near)
    // {
    //     return;
    // }

    // vec4 p0 = standardize(projection_matrix * view_matrix * pos_to_vec4(pos0));
    // vec4 p1 = standardize(projection_matrix * view_matrix * pos_to_vec4(pos1));
    // vec4 p2 = standardize(projection_matrix * view_matrix * pos_to_vec4(pos2));
    // vec4 p3 = standardize(projection_matrix * view_matrix * pos_to_vec4(pos3));
    // vec4 p4 = standardize(projection_matrix * view_matrix * pos_to_vec4(pos4));
    // vec4 p5 = standardize(projection_matrix * view_matrix * pos_to_vec4(pos5));

    array1d<vec3> points(6);
    points(0) = pos0;
    points(1) = pos1;
    points(2) = pos2;
    points(3) = pos3;
    points(4) = pos4;
    points(5) = pos5;

    array1d<vec3> curve_points = spine_curve(points.slice(), 0.05f);

    for(int i = 0; i < curve_points.size; i++)
    {
        draw_point(curve_points(i), camera, fbo, color, 0);
    }

    // draw_line(pos0, pos1, camera, fbo, color);
    // draw_line(pos1, pos2, camera, fbo, color);
    // draw_line(pos2, pos3, camera, fbo, color);
    // draw_line(pos3, pos4, camera, fbo, color);
    // draw_line(pos4, pos5, camera, fbo, color);
}

#endif
