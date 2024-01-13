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
    InitWindow(720, 720, "raylib [background]");
    SetTargetFPS(60);
    //相机
    MooCamera camera(1, 0.1f, 0.1f, 0.05f, 50.f, vec3(0.f, 0.f, 3.f), mat3(), 720, 720);
    camera.set_view(vec3(0.f, 0.f, -1.f));
    camera.update_orientation(vec3(0.f, 1.f, 0.f));
    //控制器
    Camera_Controller camera_controller;
    bind_controller(camera_controller, &camera.position, &camera.direction);
    //角色，保存原始所有信息
    Character character;
    character_load(character, "../resources/character.bin");
    //动作库
    BVH_Motion motion;
    Motion_load(motion, "../resources/long_motion.bin");
    //模拟器
    Simulator simulator;
    bind_simulator(simulator, motion, 12600, 0.05);
    //点光源
    Point_Light light(vec3(0.f, 1.f, 1.5f), vec3(0.f, 0.f, -1.f), vec3(20.f));
    //对角色数据的拷贝
    MooMesh mesh1 = make_character_rest_mesh(character);
    array1d<MooMesh> meshes(simulator.bone_masses.size);
    for(int i = 0; i < meshes.size; i++)
    {
        meshes(i) = capsule(simulator.bone_shapes(i).radius, simulator.bone_shapes(i).length, 18, 18, simulator.bone_shapes(i).radius);
    }
    //材质
    MooMaterial material1(PBR);
    material1.roughness = 0.1f;
    material1.tex = img_to_tex("../resources/skybox.png");
    material1.BRDFLut = img_to_tex("../resources/BRDFLut.png");
    material1.EavgLut = img_to_tex("../resources/EavgLut.png");
    //model存储mesh的拷贝，但由于mesh本身存储的是地址，实际上仍为mesh地址上的数据
    MooModel model0 = get_light_model(light);
    MooModel model1 = mesh_material_to_model(mesh1, material1);
    array1d<MooModel> models(simulator.bone_masses.size);
    for(int i = 0; i < models.size; i++)
    {
        models(i) = mesh_material_to_model(meshes(i), material1);
    }
    //更新信息
    updated_paramters par;
    par.light_dir = &light.light_direction;
    par.light_pos = &light.light_position;
    par.light_radiance = &light.light_radiance;
    par.shadow_map = light.shadow_map;
    //渲染器
    MooRenderer renderer;

    renderer.models.push_back(&model0);
    renderer.models.push_back(&model1);
    for(int i = 0; i < models.size; i++)
    {
        renderer.models.push_back(&models(i));
    }
    renderer.point_lights.push_back(&light);

    int key = 0;
    int frame_num = 0;
    float t = 0;
    int simulate_time = 2;

    float dt = 0.0166667f;

    int times = 900;

    array1d<vec3> target_character_bone_anim_positions(motion.nbones());
    array1d<quat> target_character_bone_anim_rotations(motion.nbones());

    array1d<vec3> curr_character_bone_anim_positions(motion.nbones());
    array1d<quat> curr_character_bone_anim_rotations(motion.nbones());

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

    cv::imshow("MOOLAB", fbo_to_img(&camera.fbo[0]));

    while (key != 27 && t < times * dt)
    {
        vec3 gamepadstick_left = gamepad_get_stick(GAMEPAD_STICK_LEFT);
        vec3 gamepadstick_right = gamepad_get_stick(GAMEPAD_STICK_RIGHT);

        BeginDrawing();
        ClearBackground(RAYWHITE);

        camera.fbo[0].set(vec3(100.f, 100.f, 200.f), 100.f);

        // batch_forward_kinematics_full(
        //     motion, 
        //     0.2f * frame_num + 12600, 
        //     target_character_bone_anim_positions, 
        //     target_character_bone_anim_rotations);

        // deform_character_anim_mesh(
        // character, 
        // target_character_bone_anim_positions, 
        // target_character_bone_anim_rotations, 
        // mesh1);

        // revise_anim_bones(simulator, target_character_bone_anim_rotations);

        // if(t >= 5 && t <= 10)
        // {
        //     simulator.bone_forces(0) = simulator.bone_forces(0) + vec3(-30.f, -30.f, 0.f);
        // }
        // else if(t >= 15 && t <= 20)
        // {
        //     simulator.bone_forces(0) = simulator.bone_forces(0) + vec3(30.f, 30.f, 0.f);
        // }
        for(int i = 0; i < simulate_time; i++)
        {
            simulator.simulate_gravity(9.8);
            simulator.simulate_damp(0.5);
            simulator.simulate_contact(0, 1, 0, 0.25, 1e-2);
            // simulator.simulate_global_control(target_character_bone_anim_rotations);
            simulator.simulate_local_control(motion.bone_local_rotations(0.5f * frame_num + 12600));
            simulator.simulate(dt / simulate_time);
        }

        deform_character_anim_bones(
            simulator, 
            curr_character_bone_anim_positions, 
            curr_character_bone_anim_rotations);

        deform_character_anim_mesh(
            character, 
            curr_character_bone_anim_positions, 
            curr_character_bone_anim_rotations, 
            mesh1);

        camera_controller.dir_gamepad_control(vec3(0.f, 0.f, 0.f), 0.1f);
        camera.update_orientation(vec3(0.f, 1.f, 0.f));
        // camera_controller.pos_gamepad_control(vec3(x(t, 0).x, 0.f, x(t, 0).z), 3);
        if(frame_num >= 0)
        {
            camera_controller.pos_gamepad_control(vec3(simulator.bone_shapes(0).pos.x, 0.f, simulator.bone_shapes(0).pos.z), camera.orientation, 3);
        }

        light.light_position = vec3(simulator.bone_shapes(0).pos.x, 3.f, simulator.bone_shapes(0).pos.z) -  
                               3.f * camera.direction;
        light.light_radiance = 20.f + 10.f * sin(deg_to_rad(0.5f * t));
        renderer.models[0]->transform = mat4(eye3(), renderer.point_lights[0]->light_position) * mat4(0.01f);

        // float delta = 0;
        for(int i = 2; i < 24; i++)
        {
            // renderer.models[i]->transform = mat4(eye3(), x(t, i - 2)) * mat4(R(t, i - 2), vec3());
            renderer.models[i]->transform = mat4(eye3(), simulator.bone_shapes(i - 2).pos) * mat4(quat_to_Rodrigues(simulator.bone_shapes(i - 2).rot), vec3());
            // delta = delta + length(x(t, i - 2) - simulator.bone_shapes(i - 2).pos);
        }
        // std::cout << t << ',' << delta << std::endl;
        
        draw(renderer.models[0], camera, &camera.fbo[0], &par, false);
        // draw(renderer.models[1], camera, &camera.fbo[0], &par, false);
        for(int i = 2; i < 24; i++)
        {
            draw(renderer.models[i], camera, &camera.fbo[0], &par, false);
        }
        draw(renderer.models[0], camera, &camera.fbo[0], &par, true);
        // draw(renderer.models[1], camera, &camera.fbo[0], &par, true);
        for(int i = 2; i < 24; i++)
        {
            draw(renderer.models[i], camera, &camera.fbo[0], &par, true);
        }

        cv::imshow("MOOLAB", fbo_to_img(&camera.fbo[0]));

        std::cout << " Time: " << t << " Frame: " << frame_num 
                  << std::endl;
        frame_num++;
        t = t + dt;
        key = cv::waitKey(1);
        EndDrawing();
    }

    CloseWindow();

    // matrix m1(195, 195);
    // matrix m2(195, 1);

    // std::ifstream infile;

    // infile.open("../resources/mattest.txt");
    // for(int i = 0; i < 195; i++)
    // {
    //     for(int j = 0; j < 195; j++)
    //     {
    //         infile >> m1(i, j);
    //     }
    // }
    // for(int i = 0; i < 195; i++)
    // {
    //     for(int j = 0; j < 1; j++)
    //     {
    //         infile >> m2(i, j);
    //     }
    // }
    // infile.close();

    // std::cout << "begin!" << std::endl;
    // for(int i = 0; i < 1; i++)
    // {
    //     matrix m3 = Elimination_method_solution_of_linear_equations(m1, m2);
    //     print(m3);
    // }
    // std::cout << "end!" << std::endl;

    // matrix m4(5, 5);
    // m4.zero();
    
    // m4(0, 0) = 1;
    // m4(0, 2) = 3;
    // m4(2, 3) = 5;
    // m4(4, 2) = 8;
    // m4(4, 4) = 3;

    // matrix m5(2, 3);
    // m5.zero();
    // m5(0, 0) = 1.43214;
    // m5(0, 1) = -2.4123;
    // m5(0, 2) = -1.2134;
    // m5(1, 2) = -3.43214;

    // vec3 m6(1.f, 2.f, 5.f);

    // mat3 m7 = eye3();

    // set_matrix(m4, 1, 2, 2, 4, m5);
    // set_matrix(m4, 1, 3, 1, 1, m6);
    // set_matrix(m4, 2, 4, 2, 4, m7);

    // print(m4);

    // m4 = m4 * m4 - m4 + m4 * transpose_matrix(m4) * m4;

    // print(m4);

    // vec3 m8;

    // set_vec(m8, m4, 2, 4, 2, 2);

    // print(m8);

    // print(cross(m6, m8));

    // print(Rodrigues(24, cross(m6, m8)));

    // print(quat(24, cross(m6, m8)));

    // print(quat_to_Rodrigues(quat(24, cross(m6, m8))));

    // print(Rodrigues(rad_to_deg(length(m6) * 0.1), m6) - quat_to_Rodrigues(avel_to_quat(m6, 0.1)));

    return 0;
}
