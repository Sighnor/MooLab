#ifndef MOOLAB_STRATEGY
#define MOOLAB_STRATEGY

#include "array.hpp"
#include "global.hpp"
#include "vec.hpp"

struct Code
{
    int num;
    vec2 coord;
};

struct Graph
{
    array1d<Code> codes;

    int size() const {return codes.size;}

    float &coord_x()
    {
        assert(codes.size == 16);
        return codes(15).coord.x;
    }

    float &coord_y()
    {
        assert(codes.size == 16);
        return codes(15).coord.y;
    }

    int get_code_num(int id)
    {
        assert(id >= 0 && id < codes.size);
        return codes(id).num;
    }
    vec2 get_code_world_coord(int id)
    {
        assert(id >= 0 && id < codes.size);
        return vec2((4 - codes(id).coord.y) * 720 / 5, (1 + codes(id).coord.x) * 720 / 5);
    }

    float get_Astar_cost()
    {
        float cost = 0;
        int id = 0;
        for(int i = 0; i < 4; i++)
        {
            for(int j = 0; j < 4 && id < 15; j++)
            {
                cost = cost + abs(codes(id).coord.x - i) + abs(codes(id).coord.y - j);
                id++;
            }
        }
        return cost;
    }
};

Graph init_graph()
{
    Graph graph;
    graph.codes.resize(16);
    for(int i = 0; i < 4; i++)
    {
        for(int j = 0; j < 4; j++)
        {
            int id = i * 4 + j;
            graph.codes(id).num = id + 1;
            graph.codes(id).coord = vec2(i, j);
        }
    }
    return graph;
}

Graph update_graph(Graph &graph, int dir)
{
    Graph updated_graph;
    updated_graph.codes = graph.codes;
    if (dir == 0)
    {
        updated_graph.coord_x() = updated_graph.coord_x() - 1;
    }
    else if (dir == 1)
    {
        updated_graph.coord_x() = updated_graph.coord_x() + 1;
    }
    else if (dir == 2)
    {
        updated_graph.coord_y() = updated_graph.coord_y() - 1;
    }
    else if (dir == 3)
    {
        updated_graph.coord_y() = updated_graph.coord_y() + 1;
    }

    for(int i = 0; i < updated_graph.size() - 1; i++)
    {
        if(length(graph.codes(i).coord - updated_graph.codes(15).coord) < 1e-4f)
        {
            updated_graph.codes(i).coord = graph.codes(15).coord;
        }
    }
    return updated_graph;
}

struct Node
{
    bool closed;
    int depth;
    int last_dir;
    int parent_id;
    Graph graph;

    float get_total_cost()
    {
        return depth + graph.get_Astar_cost();
    }

    float get_Astar_cost()
    {
        return graph.get_Astar_cost();
    }

    bool can_extand(int dir)
    {
        if(dir == 0 && last_dir!= 1 && graph.coord_x() > 0)
        {
            return true;
        }
        else if(dir == 1 && last_dir!= 0 && graph.coord_x() < 3)
        {
            return true;
        }
        else if(dir == 2 && last_dir!= 3 && graph.coord_y() > 0)
        {
            return true;
        }
        else if(dir == 3 && last_dir!= 2 && graph.coord_y() < 3)
        {
            return true;
        }
        return false;
    }
};

Node init_node(Graph &graph)
{
    Node node;
    node.closed = false;
    node.depth = 0;
    node.last_dir = -1;
    node.parent_id = -1;
    node.graph = graph;

    return node;
}

Node extand_node(Node &node, int id, int dir)
{
    Node extanded_node;
    extanded_node.closed = false;
    extanded_node.depth = node.depth + 1;
    extanded_node.last_dir = dir;
    extanded_node.parent_id = id;
    extanded_node.graph = update_graph(node.graph, dir);

    return extanded_node;
}

std::vector<Node> Astar(Graph &graph, int &end_id)
{
    std::vector<Node> nodes;

    Node root_node = init_node(graph);

    nodes.push_back(root_node);

    int t = 0;
    int open_id = 0;
    float least_Astar_cost = 10000;
    float least_total_cost = 10000;

    while (least_Astar_cost != 0 && t < 10000)
    {
        for(int i = 0; i < 4; i++)
        {
            if(nodes[open_id].can_extand(i))
            {
                Node extanded_node = extand_node(nodes[open_id], open_id, i);
                nodes.push_back(extanded_node);
            }
        }
        nodes[open_id].closed = true;

        least_total_cost = 10000;
        for(int i = 0; i < nodes.size(); i++)
        {
            if(nodes[i].get_Astar_cost() < least_Astar_cost)
            {
                least_Astar_cost = nodes[i].get_Astar_cost();
                end_id = i;
            }
            if(nodes[i].closed == false && nodes[i].get_total_cost() < least_total_cost)
            {
                least_total_cost = nodes[i].get_total_cost();
                open_id = i;
            }
        }

        t++;
        std::cout << t << ", " << least_Astar_cost << std::endl;
    }

    return nodes;
}

#endif