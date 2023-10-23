#ifndef MOOLAB_UI
#define MOOLAB_UI

#include "texture.hpp"
#include "vec.hpp"

struct Moochar
{
    char name;
    std::vector<vec2> points;
};

struct MooUI
{
    int size;
    std::map<char, int> distances;
    std::vector<Moochar> moochars;
};

void UI_load_char(MooUI &ui, char c, const char* img_name)
{
    Moochar moochar;
    moochar.name = c;

    texture tex = img_to_tex(img_name);

    int max_i = 0;
    int min_i = tex.rows;
    int max_j = 0;
    int min_j = tex.cols;

    for(int i = 0; i < tex.rows; i++)
    {
        for(int j = 0; j < tex.cols; j++)
        {
            if(length(tex.getcolor(vec2(float(i) / tex.rows, float(j) / tex.cols))) < 400.f)
            {
                if(i >= max_i)
                {
                    max_i = i;
                }
                if(i <= min_i)
                {
                    min_i = i;
                }
                if(j >= max_j)
                {
                    max_j = j;
                }
                if(j <= min_j)
                {
                    min_j = j;
                }

                moochar.points.push_back(vec2(i, j));
            }
        }
    }

    vec2 center_point = vec2((min_i + max_i) / 2, (min_j + max_j) / 2);

    for(auto &p : moochar.points)
    {
        p = p - center_point;
    }

    ui.distances[c] = ui.size;
    ui.moochars.push_back(moochar);
    ui.size = ui.size + 1;

    // std::cout << c << ',' << moochar.points.size() << std::endl;
}

void UI_draw_char(char c, vec2 point, MooUI &ui, FBO *fbo, vec3 color)
{  
    int id;
    id = ui.distances[c];
    if(c >= 48 && c <= 57)
    {
        id = c - 48;
    }
    else if(c >= 65 && c <= 90)
    {
        id = c - 65 + 10;
    }
    else if(c >= 97 && c <= 122)
    {
        id = c - 97 + 10;
    }
    assert(id >= 0 && id < ui.size);
    for(auto p : ui.moochars[id].points)
    {
        int i = point.x + p.x;
        int j = point.y + p.y;

        if(i >= 0 && i <= fbo->rows - 1 && j >= 0 && j <= fbo->cols - 1)
        {
            fbo->getcolor(i, j) = clampv(color);
        }
    }
}

void UI_draw_num(int num, vec2 point, MooUI &ui, FBO *fbo, vec3 color)
{  
    int id;
    id = num;
    if(id >= 0 && id < 10)
    {
        for(auto p : ui.moochars[id].points)
        {
            int i = point.x + p.x;
            int j = point.y + p.y;

            if(i >= 0 && i <= fbo->rows - 1 && j >= 0 && j <= fbo->cols - 1)
            {
                fbo->getcolor(i, j) = clampv(color);
            }
        }
    }
    if(id >= 10 && id <= 100)
    {
        int id0 = id % 10;
        int id1 = id / 10;

        vec2 point0 = vec2(point.x, point.y + 25);
        vec2 point1 = vec2(point.x, point.y - 25);

        for(auto p : ui.moochars[id0].points)
        {
            int i = point0.x + p.x;
            int j = point0.y + p.y;

            if(i >= 0 && i <= fbo->rows - 1 && j >= 0 && j <= fbo->cols - 1)
            {
                fbo->getcolor(i, j) = clampv(color);
            }
        }

        for(auto p : ui.moochars[id1].points)
        {
            int i = point1.x + p.x;
            int j = point1.y + p.y;

            if(i >= 0 && i <= fbo->rows - 1 && j >= 0 && j <= fbo->cols - 1)
            {
                fbo->getcolor(i, j) = clampv(color);
            }
        }
    }
}

#endif