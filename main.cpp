extern "C"
{
#include "raylib.h"
#include "raymath.h"
}

#include "motion.hpp"
#include "character.hpp"
#include "curve.hpp"
#include "database.hpp"
#include "ik.hpp"
#include "motion.hpp"
#include "nnet.hpp"
#include "renderer.hpp"
#include "robot.hpp"
#include "simulator.hpp"

enum
{
    GAMEPAD_PLAYER = 0,
};

enum
{
    GAMEPAD_STICK_LEFT,
    GAMEPAD_STICK_RIGHT,
};

vec3 gamepad_get_stick(int stick, const float deadzone = 0.2f)
{
    float gamepadx = GetGamepadAxisMovement(GAMEPAD_PLAYER, stick == GAMEPAD_STICK_LEFT ? GAMEPAD_AXIS_LEFT_X : GAMEPAD_AXIS_RIGHT_X);
    float gamepady = GetGamepadAxisMovement(GAMEPAD_PLAYER, stick == GAMEPAD_STICK_LEFT ? GAMEPAD_AXIS_LEFT_Y : GAMEPAD_AXIS_RIGHT_Y);
    float gamepadmag = sqrtf(gamepadx * gamepadx + gamepady * gamepady);
    
    if (gamepadmag > deadzone)
    {
        float gamepaddirx = gamepadx / gamepadmag;
        float gamepaddiry = gamepady / gamepadmag;
        float gamepadclippedmag = gamepadmag > 1.0f ? 1.0f : gamepadmag * gamepadmag;
        gamepadx = gamepaddirx * gamepadclippedmag;
        gamepady = gamepaddiry * gamepadclippedmag;
    }
    else
    {
        gamepadx = 0.0f;
        gamepady = 0.0f;
    }
    
    return vec3(gamepadx, 0.0f, gamepady);
}

int main(int argc, char** argv)
{
    // InitWindow(720, 720, "raylib [background]");
    // SetTargetFPS(60);
    // //相机
    // MooCamera camera(1, 0.1f, 0.1f, 0.05f, 50.f, vec3(0.f, 0.f, 3.f), mat3(), 720, 720);
    // camera.set_view(vec3(0.f, 0.f, -1.f));
    // camera.update_orientation(vec3(0.f, 1.f, 0.f));
    // //控制器
    // Camera_Controller camera_controller;
    // bind_controller(camera_controller, &camera.position, &camera.direction);
    // //点光源
    // Point_Light light(vec3(0.f, 1.f, 1.5f), vec3(0.f, 0.f, -1.f), vec3(20.f));
    // //对角色数据的拷贝
    // MooMesh mesh1 = capsule(0.5f * 0.0738f, 0.0738f, 18, 18,  0.5f * 0.0738f);
    // MooMesh mesh2 = capsule(0.5f * 0.1f, 0.435f, 18, 18,  0.5f * 0.1f);
    // MooMesh mesh3 = capsule(0.5f * 0.1f, 0.4237f, 18, 18,  0.5f * 0.1f);
    // MooMesh mesh4 = capsule(0.5f * 0.1f, 0.1730f, 18, 18,  0.5f * 0.1f);
    // MooMesh mesh5 = capsule(0.5f * 0.05f, 0.05f, 18, 18,  0.5f * 0.05f);
    // MooMesh mesh6 = capsule(0.5f * 0.12f, 0.1259f, 18, 18,  0.5f * 0.12f);
    // MooMesh mesh7 = capsule(0.5f * 0.12f, 0.1234f, 18, 18,  0.5f * 0.12f);
    // MooMesh mesh8 = capsule(0.5f * 0.12f, 0.2583f, 18, 18,  0.5f * 0.12f);
    // MooMesh mesh9 = capsule(0.5f * 0.10f, 0.1177f, 18, 18,  0.5f * 0.10f);
    // MooMesh mesh10 = capsule(0.5f * 0.15f, 0.25f, 18, 18,  0.5f * 0.15f);
    // MooMesh mesh11 = capsule(0.5f * 0.1128f, 0.2583f, 18, 18,  0.5f * 0.1128f);
    // MooMesh mesh12 = capsule(0.5f * 0.1f, 0.33f, 18, 18,  0.5f * 0.1f);
    // MooMesh mesh13 = capsule(0.5f * 0.1f, 0.252f, 18, 18,  0.5f * 0.1f);
    // MooMesh mesh14 = capsule(0.5f * 0.1f, 0.2f, 18, 18,  0.5f * 0.1f);
    // //材质
    // MooMaterial material1(TEXTURE);
    // material1.roughness = 0.1f;
    // material1.tex = img_to_tex("../resources/skybox.png");
    // material1.BRDFLut = img_to_tex("../resources/BRDFLut.png");
    // material1.EavgLut = img_to_tex("../resources/EavgLut.png");
    // //model存储mesh的拷贝，但由于mesh本身存储的是地址，实际上仍为mesh地址上的数据
    // MooModel model0 = get_light_model(light);
    // MooModel model1 = mesh_material_to_model(mesh1, material1);
    // MooModel model2 = mesh_material_to_model(mesh2, material1);
    // MooModel model3 = mesh_material_to_model(mesh3, material1);
    // MooModel model4 = mesh_material_to_model(mesh4, material1);
    // MooModel model5 = mesh_material_to_model(mesh5, material1);
    // MooModel model6 = mesh_material_to_model(mesh2, material1);
    // MooModel model7 = mesh_material_to_model(mesh3, material1);
    // MooModel model8 = mesh_material_to_model(mesh4, material1);
    // MooModel model9 = mesh_material_to_model(mesh5, material1);
    // MooModel model10 = mesh_material_to_model(mesh6, material1);
    // MooModel model11 = mesh_material_to_model(mesh7, material1);
    // MooModel model12 = mesh_material_to_model(mesh8, material1);
    // MooModel model13 = mesh_material_to_model(mesh9, material1);
    // MooModel model14 = mesh_material_to_model(mesh10, material1);
    // MooModel model15 = mesh_material_to_model(mesh11, material1);
    // MooModel model16 = mesh_material_to_model(mesh12, material1);
    // MooModel model17 = mesh_material_to_model(mesh13, material1);
    // MooModel model18 = mesh_material_to_model(mesh14, material1);
    // MooModel model19 = mesh_material_to_model(mesh11, material1);
    // MooModel model20 = mesh_material_to_model(mesh12, material1);
    // MooModel model21 = mesh_material_to_model(mesh13, material1);
    // MooModel model22 = mesh_material_to_model(mesh14, material1);
    // //更新信息
    // updated_paramters par;
    // par.light_dir = &light.light_direction;
    // par.light_pos = &light.light_position;
    // par.light_radiance = &light.light_radiance;
    // par.shadow_map = light.shadow_map;
    // //渲染器
    // MooRenderer renderer;

    // renderer.models.push_back(&model0);
    // renderer.models.push_back(&model1);
    // renderer.models.push_back(&model2);
    // renderer.models.push_back(&model3);
    // renderer.models.push_back(&model4);
    // renderer.models.push_back(&model5);
    // renderer.models.push_back(&model6);
    // renderer.models.push_back(&model7);
    // renderer.models.push_back(&model8);
    // renderer.models.push_back(&model9);
    // renderer.models.push_back(&model10);
    // renderer.models.push_back(&model11);
    // renderer.models.push_back(&model12);
    // renderer.models.push_back(&model13);
    // renderer.models.push_back(&model14);
    // renderer.models.push_back(&model15);
    // renderer.models.push_back(&model16);
    // renderer.models.push_back(&model17);
    // renderer.models.push_back(&model18);
    // renderer.models.push_back(&model19);
    // renderer.models.push_back(&model20);
    // renderer.models.push_back(&model21);
    // renderer.models.push_back(&model22);
    // renderer.point_lights.push_back(&light);

    // int key = 0;
    // int frame_num = 0;
    // int t = 0;

    // float dt = 0.0166667f;

    // int times = 3000;

    // array2d<vec3> x(times, 22);
    // array2d<mat3> R(times, 22);

    // std::ifstream infile;

    // infile.open("../resources/x.txt");
    // for(int i = 0; i < times; i++)
    // {
    //     for(int j = 0; j < 22; j++)
    //     {
    //         infile >> x(i, j).x;
    //     }
    //     for(int j = 0; j < 22; j++)
    //     {
    //         infile >> x(i, j).y;
    //     }
    //     for(int j = 0; j < 22; j++)
    //     {
    //         infile >> x(i, j).z;
    //     }
    // }
    // infile.close();

    // infile.open("../resources/R.txt");
    // for(int i = 0; i < times; i++)
    // {
    //     for(int j = 0; j < 22; j++)
    //     {
    //         infile >> R(i, j).X.x;
    //         infile >> R(i, j).Y.x;
    //         infile >> R(i, j).Z.x;
    //     }
    //     for(int j = 0; j < 22; j++)
    //     {
    //         infile >> R(i, j).X.y;
    //         infile >> R(i, j).Y.y;
    //         infile >> R(i, j).Z.y;
    //     }
    //     for(int j = 0; j < 22; j++)
    //     {
    //         infile >> R(i, j).X.z;
    //         infile >> R(i, j).Y.z;
    //         infile >> R(i, j).Z.z;
    //     }
    // }
    // infile.close();

    // cv::imshow("MOOLAB", fbo_to_img(&camera.fbo[0]));

    // while (key != 27 && t < times)
    // {
    //     vec3 gamepadstick_left = gamepad_get_stick(GAMEPAD_STICK_LEFT);
    //     vec3 gamepadstick_right = gamepad_get_stick(GAMEPAD_STICK_RIGHT);

    //     BeginDrawing();
    //     ClearBackground(RAYWHITE);

    //     camera.fbo[0].set(vec3(100.f, 100.f, 200.f), 100.f);

    //     // camera.position = camera.position + 5.f * (quat(rad_to_deg(camera_controller.ang.x), vec3(0.f, 1.f, 0.f)) * (-gamepadstick_left)) * dt;
    //     // camera_controller.dir_gamepad_control(gamepadstick_right, 0.1f);
    //     camera_controller.dir_gamepad_control(vec3(0.05f, 0.f, 0.f), 0.1f);
    //     camera.update_orientation(vec3(0.f, 1.f, 0.f));
    //     camera.position = vec3(0.f, 0.f, 0.f + x(t, 0).z) + quat(rad_to_deg(camera_controller.ang.x), vec3(0.f, 1.f, 0.f)) * vec3(0.f, 0.f, -3.f);
    //     // camera.position = x(t, 0) + quat(rad_to_deg(camera_controller.ang.x), vec3(0.f, 1.f, 0.f)) * vec3(0.f, 0.f, -3.f);
    //     // camera_controller.pos_pid_control(1.f, 0.05f, 0.1f, 500.f, dt, vec3(x(t, 0).x, x(t, 0).y, x(t, 0).z + 3.f));

    //     renderer.models[0]->transform = mat4(eye3(), renderer.point_lights[0]->light_position) * mat4(0.01f);
    //     for(int i = 1; i < 23; i++)
    //     {
    //         if((i >= 2 && i <= 5) || (i >= 15 && i <= 18))
    //         {
    //             renderer.models[i]->transform = mat4(Rodrigues(0, vec3(0.f, 1.f, 0.f)), vec3()) *  mat4(eye3(), x(t, i - 1)) * mat4(R(t, i - 1), vec3());
    //         }
    //         else if((i >= 6 && i <= 9) || (i >= 19 && i <= 22))
    //         {
    //             renderer.models[i]->transform = mat4(Rodrigues(0, vec3(0.f, 1.f, 0.f)), vec3()) *  mat4(eye3(), x(t, i - 1)) * mat4(R(t, i - 1), vec3());
    //         }
    //         else
    //         {
    //             renderer.models[i]->transform = mat4(Rodrigues(0, vec3(0.f, 1.f, 0.f)), vec3()) *  mat4(eye3(), x(t, i - 1)) * mat4(R(t, i - 1), vec3());
    //         }
    //     }
        
    //     draw(renderer.models[0], camera, &camera.fbo[0], &par, false);
    //     for(int i = 1; i < 23; i++)
    //     {
    //         draw(renderer.models[i], camera, &camera.fbo[0], &par, false);
    //     }
    //     draw(renderer.models[0], camera, &camera.fbo[0], &par, true);
    //     for(int i = 1; i < 23; i++)
    //     {
    //         draw(renderer.models[i], camera, &camera.fbo[0], &par, true);
    //     }

    //     cv::imshow("MOOLAB", fbo_to_img(&camera.fbo[0]));

    //     std::cout << " Time: " << t << " Frame: " << frame_num 
    //               << std::endl;
    //     frame_num++;
    //     t++;
    //     key = cv::waitKey(1);
    //     EndDrawing();
    // }

    // CloseWindow();

    matrix m1(195, 195);
    matrix m2(195, 1);

    std::ifstream infile;

    infile.open("../resources/mattest.txt");
    for(int i = 0; i < 195; i++)
    {
        for(int j = 0; j < 195; j++)
        {
            infile >> m1(i, j);
        }
    }
    for(int i = 0; i < 195; i++)
    {
        for(int j = 0; j < 1; j++)
        {
            infile >> m2(i, j);
        }
    }
    infile.close();

    std::cout << "begin!" << std::endl;
    for(int i = 0; i < 1000; i++)
    {
        matrix m3 = Elimination_method_solution_of_linear_equations(m1, m2);
    }
    std::cout << "end!" << std::endl;

    return 0;
}
