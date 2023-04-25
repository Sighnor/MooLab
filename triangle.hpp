#ifndef ENGINE_TRIANGLE
#define ENGINE_TRIANGLE

#include "vec.hpp"

class Triangle
{
    public:
        Triangle();

        vec3 pos[3];
        vec3 normal[3];
        vec2 texcoords[3];
        vec3 color[3];

        void set_vertex(int id, vec3 ver);
        void set_normal(int id, vec3 n);
        void set_texcoords(int id, vec2 uv);
        void set_color(int id, vec3 col);
};

Triangle::Triangle()
{
    pos[0] = vec3();
    pos[1] = vec3();
    pos[2] = vec3();

    color[0] = vec3();
    color[1] = vec3();
    color[2] = vec3();

    texcoords[0] = vec2();
    texcoords[1] = vec2();
    texcoords[2] = vec2();

    normal[0] = vec3();
    normal[1] = vec3();
    normal[2] = vec3(); 
}

void Triangle::set_vertex(int id, vec3 ver)
{
    pos[id] = ver;
}

void Triangle::set_normal(int id, vec3 n)
{
    normal[id] = n;
}

void Triangle::set_texcoords(int id, vec2 uv)
{
    texcoords[id] = uv;
}

void Triangle::set_color(int id, vec3 col)
{
    color[id] = col;
}

#endif