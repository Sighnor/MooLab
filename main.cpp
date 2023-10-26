#include "strategy.hpp"
#include "ui.hpp"


int main(int argc, char** argv)
{
    //codes
    Graph graph = init_graph();
    //UI
    MooUI ui;
    UI_load_char(ui, '0', "../resources/UI/0.png");
    UI_load_char(ui, '1', "../resources/UI/1.png");
    UI_load_char(ui, '2', "../resources/UI/2.png");
    UI_load_char(ui, '3', "../resources/UI/3.png");
    UI_load_char(ui, '4', "../resources/UI/4.png");
    UI_load_char(ui, '5', "../resources/UI/5.png");
    UI_load_char(ui, '6', "../resources/UI/6.png");
    UI_load_char(ui, '7', "../resources/UI/7.png");
    UI_load_char(ui, '8', "../resources/UI/8.png");
    UI_load_char(ui, '9', "../resources/UI/9.png");

    FBO output(720, 720);

    int key = 0;
    int t = 0;
    int velid;
    vec2 vel;

    while (key != 27)
    {
        output.set(vec3(255.f, 255.f, 255.f), 100.f);

        for(int i = 0; i < graph.size() - 1; i++)
        {
            UI_draw_num(graph.get_code_num(i), graph.get_code_world_coord(i), ui, &output, vec3(0.f, 0.f, 0.f));
        }

        cv::imshow("MOOLAB", fbo_to_img(&output));

        if (t % 1 == 0) 
        {
            vel = vec2(0.f);

            vec2 temp_coord = graph.codes(15).coord;
            if (key == 97 && graph.coord_x() > 0)
            {

                graph = update_graph(graph, 0);
                // graph.codes(15).coord.x = clampf(graph.codes(15).coord.x - 1, 0, 3);
            }
            else if (key == 100 && graph.coord_x() < 3)
            {
                graph = update_graph(graph, 1);
                // graph.codes(15).coord.x = clampf(graph.codes(15).coord.x + 1, 0, 3);
            }
            else if (key == 115 && graph.coord_y() > 0)
            {
                graph = update_graph(graph, 2);
                // graph.codes(15).coord.y = clampf(graph.codes(15).coord.y - 1, 0, 3);
            }
            else if (key == 119 && graph.coord_y() < 3)
            {
                graph = update_graph(graph, 3);
                // graph.codes(15).coord.y = clampf(graph.codes(15).coord.y + 1, 0, 3);
            }

            std::cout << graph.coord_x() << ", " << graph.coord_y() << std::endl;

            // for(int i = 0; i < graph.size() - 1; i++)
            // {
            //     if(length(graph.codes(i).coord - graph.codes(15).coord) < 1e-4f)
            //     {
            //         velid = i;
            //         vel = (temp_coord - graph.codes(15).coord) / 29.f;
            //         graph.codes(15).coord = graph.codes(i).coord;
            //         break;
            //     }
            // }
        }
        // else
        // {
        //     graph.codes(velid).coord = graph.codes(velid).coord + vel;
        // }

        // std::cout << " Time: " << t << std::endl;
        t++;
        key = cv::waitKey(10);
    }

    // std::vector<Node> nodes;

    // Node root_node = init_node(graph);

    // nodes.push_back(root_node);

    // int open_id = 0;
    // int end_id = 0;
    // float least_total_cost = 10000;
    // float least_Astar_cost = 10000;

    // while (least_Astar_cost != 0)
    // {
    //     least_total_cost = 10000;
    //     for(int i = 0; i < nodes.size(); i++)
    //     {
    //         if(nodes[i].get_Astar_cost() < least_Astar_cost)
    //         {
    //             least_Astar_cost = nodes[i].get_Astar_cost();
    //             end_id = i;
    //         }
    //         if(nodes[i].closed == false && nodes[i].get_total_cost() < least_total_cost)
    //         {
    //             least_total_cost = nodes[i].get_total_cost();
    //             open_id = i;
    //         }
    //     }

    //     for(int i = 0; i < 4; i++)
    //     {
    //         if(nodes[open_id].can_extand(i))
    //         {
    //             Node extanded_node = extand_node(nodes[open_id], open_id, i);
    //             nodes.push_back(extanded_node);
    //         }
    //     }
    //     nodes[open_id].closed = true;

    //     std::cout << open_id << ',' << least_Astar_cost << std::endl;
    // }

    int end_id = 0;
    std::vector<Node> nodes = Astar(graph, end_id);

    int id = end_id;
    std::vector<int> ids;
    while(id != -1)
    {
        std::cout << nodes[id].depth << std::endl;
        ids.push_back(id);
        id = nodes[id].parent_id;
    }

    for(int j = ids.size() - 1; j >= 0; j--)
    {
        output.set(vec3(255.f, 255.f, 255.f), 100.f);

        for(int i = 0; i < nodes[ids[j]].graph.size() - 1; i++)
        {
            UI_draw_num(nodes[ids[j]].graph.get_code_num(i), nodes[ids[j]].graph.get_code_world_coord(i), ui, &output, vec3(0.f, 0.f, 0.f));
        }

        cv::imshow("MOOLAB", fbo_to_img(&output));
        key = cv::waitKey(10000);
    }

    return 0;
}
