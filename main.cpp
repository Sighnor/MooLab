extern "C"
{
#include "raylib.h"
#include "raymath.h"
}

#include "character.hpp"
#include "curve.hpp"
#include "database.hpp"
#include "ik.hpp"
#include "motion.hpp"
#include "nnet.hpp"
#include "renderer.hpp"
#include "robot.hpp"
#include "simulator.hpp"
#include "ui.hpp"

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
    MooCamera camera(1, 0.1f, 0.1f, 0.05f, 50.f, vec3(0.f, 1.3f, -6.f), mat3(), 720, 720);
    camera.set_view(vec3(0.f, 0.f, 1.f));
    camera.update_orientation(vec3(0.f, 1.f, 0.f));
    //控制器
    Camera_Controller camera_controller;
    bind_controller(camera_controller, &camera.position, &camera.direction);
    Character_Controller character_controller;
    character_controller.rotation0 = quat(1.f, 0.f, 0.f, 0.f);
    //角色，保存原始所有信息
    Character character;
    character_load(character, "../resources/character.bin");
    //数据库
    Database db;
    Database_load(db, "../resources/database.bin");
    //动作库
    int N = 50;
    int dN = 10;
    BVH_Motion motion, motion1, motion2, motion3, motion4, motion5, motion6, active_motion, search_motion;
    Motion_load(motion, "../resources/long_motion.bin");
    motion1 = motion_sub_sequence(motion, 0, 281);
    motion2 = motion_sub_sequence(motion, 281, 562);
    motion3 = motion_sub_sequence(motion, 562, 13153);
    motion4 = motion_sub_sequence(motion, 13153, 25744);
    motion5 = motion_sub_sequence(motion, 25744, 39622);
    motion6 = motion_sub_sequence(motion, 39622, 53500);
    motion = concatenate_two_motions(motion1, motion2, 281, 30);
    motion = concatenate_two_motions(motion, motion3, 562, 30);
    motion = concatenate_two_motions(motion, motion4, 13153, 30);
    motion = concatenate_two_motions(motion, motion5, 25744, 30);
    motion = concatenate_two_motions(motion, motion6, 39622, 30);
    active_motion = motion_sub_sequence(motion, 0, 0 + 1.5 * N);
    search_motion = motion_sub_sequence(motion, 0, 0 + 1.5 * N);
    //模拟器
    Simulator simulator;
    bind_simulator(simulator, motion, 0, 0.05);
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
    MooMaterial material0(PBR);
    material0.roughness = 0.1f;
    material0.tex = img_to_tex("../resources/skybox.png");
    material0.BRDFLut = img_to_tex("../resources/BRDFLut.png");
    material0.EavgLut = img_to_tex("../resources/EavgLut.png");
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
    int t = 0;
    int simulate_time = 1;

    float dt = 0.0166667f;

    int times = 900;

    array1d<vec3> target_character_bone_anim_positions(motion.nbones());
    array1d<quat> target_character_bone_anim_rotations(motion.nbones());

    array1d<vec3> curr_character_bone_anim_positions(motion.nbones());
    array1d<quat> curr_character_bone_anim_rotations(motion.nbones());

    array1d<vec3> last_character_bone_anim_positions(motion.nbones());
    array1d<quat> last_character_bone_anim_rotations(motion.nbones());

    array1d<vec3> next_character_bone_anim_positions(motion.nbones());
    array1d<quat> next_character_bone_anim_rotations(motion.nbones());
    
    array1d<int> frames(motion.nframes() - N);
    array1d<int> contacts(motion.nframes() - N);
    // frame 0.0 * N 0.2 * N 0.4 * N 0.6 * N 0.8 * N 1.0 * N
    array2d<vec3> vels(motion.nframes() - N, 6);
    array2d<vec3> avels(motion.nframes() - N, 6);
    // ltoe rtoe lvel rvel
    array2d<vec3> features(motion.nframes() - N, 4);
    vec3 vel;
    vec3 avel;
    vec3 ltoe;
    vec3 rtoe;
    vec3 lvel;
    vec3 rvel;
    float contact = 0;
    vec3 contact_point;

    int weight_vel = 8;
    int weight_avel = 3;
    int weight_toe = 5;

    batch_forward_kinematics_full(
        motion, 
        0, 
        curr_character_bone_anim_positions, 
        curr_character_bone_anim_rotations);
    contact = 0;
    contact_point = curr_character_bone_anim_positions(Bone_LeftToe);

    // std::cout << db.frames.size << ',' << db.velocities.size << ',' << db.features.rows << ',' << db.features.cols << std::endl;

    cv::imshow("MOOLAB", fbo_to_img(&camera.fbo[0]));
    cv::createTrackbar("vel", "MOOLAB", &weight_vel, 50);
    cv::createTrackbar("avel", "MOOLAB", &weight_avel, 50);
    cv::createTrackbar("toe", "MOOLAB", &weight_toe, 50);

    // active_motion = motion;
    active_motion = motion_sub_sequence(motion, 12600, 13600);

    while (key != 27 && t < motion.nframes() - N)
    {
        vec3 gamepadstick_left = gamepad_get_stick(GAMEPAD_STICK_LEFT);
        vec3 gamepadstick_right = gamepad_get_stick(GAMEPAD_STICK_RIGHT);

        vec3 global_gamepad_vel = 5.f * (quat(rad_to_deg(camera_controller.ang.x), vec3(0.f, 1.f, 0.f)) * (-gamepadstick_left));
        vec3 global_gamepad_avel = cross(normalize(curr_character_bone_anim_rotations(Bone_Entity) * vec3(0.f, 0.f, 1.f)), normalize(global_gamepad_vel)) / (N * 0.0166667f);

        vec3 local_gamepad_vel = (inv_quat(curr_character_bone_anim_rotations(Bone_Entity)) * global_gamepad_vel);
        vec3 local_gamepad_avel = global_gamepad_avel;

        BeginDrawing();
        ClearBackground(RAYWHITE);

        camera.fbo[0].set(vec3(100.f, 100.f, 200.f), 100.f);

        // if(frame_num == N)
        // {
            // vel = (2.5f + 2.5f * sin(deg_to_rad(0.1f * t))) * vec3(0.f, 0.f, 1.f);

            // avel = vec3(0.f, 0.f, 0.f);

            // std::cout << 2.5f + 2.5f * sin(deg_to_rad(0.1f * t)) << ',' << vel.x << ',' << vel.y << ',' << vel.z << std::endl;
        // }
        
        // if(frame_num == N)
        // {
        //     frame_num = 0;

        //     // character_controller.vel5 = vel;
        //     // character_controller.avel5 = avel;
        //     character_controller.update(local_gamepad_vel, local_gamepad_avel, curr_character_bone_anim_rotations(Bone_Entity));
        //     int best_frame = Database_search(db, character_controller, contact, ltoe, rtoe, N, dN, weight_vel, weight_avel, weight_toe);
        //     search_motion = motion_sub_sequence(
        //                             motion, 
        //                             best_frame, 
        //                             best_frame + 1.5 * N);
        //     active_motion = motion_sub_sequence(
        //                         concatenate_two_motions(
        //                             active_motion, 
        //                             search_motion, 
        //                             N, 
        //                             0.5 * N), 
        //                         N,
        //                         N + 1.5 * N);
        // }

        // if(frame_num >= 50 && length(local_gamepad_vel) > 0.1f && (length(local_gamepad_vel - vel) / length(local_gamepad_vel) > 0.1f) || frame_num == N)
        // {
        //     character_controller.update(local_gamepad_vel, local_gamepad_avel, curr_character_bone_anim_rotations(Bone_Entity));
        //     int best_frame = Database_search(db, character_controller, contact, ltoe, rtoe, N, dN, weight_vel, weight_avel, weight_toe);
        //     search_motion = motion_sub_sequence(
        //                             motion, 
        //                             best_frame, 
        //                             best_frame + 1.5 * N);
        //     active_motion = motion_sub_sequence(
        //                         concatenate_two_motions(
        //                             active_motion, 
        //                             search_motion, 
        //                             frame_num, 
        //                             0.5 * N), 
        //                         frame_num,
        //                         frame_num + 1.5 * N);
            
        //     frame_num = 0;
        // }

        last_character_bone_anim_positions = curr_character_bone_anim_positions;
        last_character_bone_anim_rotations = curr_character_bone_anim_rotations;

        if(t <= -1)
        {
            batch_forward_kinematics_full(
            active_motion, 
            frame_num, 
            curr_character_bone_anim_positions, 
            curr_character_bone_anim_rotations);

            if(t == 1000)
            {
                // bind_simulator(simulator, active_motion, frame_num, 0.05);
                update_simulator(simulator, 
                    last_character_bone_anim_positions, 
                    last_character_bone_anim_rotations, 
                    curr_character_bone_anim_positions, 
                    curr_character_bone_anim_rotations, 
                    0.05, 
                    dt);
            }
        }
        else
        {
            for(int i = 0; i < simulate_time; i++)
            {
                simulator.bone_forces(0) = simulator.bone_forces(0) + vec3(100.f, 0.f, 0.f);
                simulator.simulate_gravity(9.8);
                simulator.simulate_damp(0.25);
                simulator.simulate_contact(0, 1, 0, 0.05, 1e-2);
                // simulator.simulate_global_control(active_motion.bone_local_rotations(0.5f * frame_num));
                // simulator.simulate_local_control(active_motion.bone_local_rotations(0.5f * frame_num));
                simulator.simulate(dt / simulate_time);
            }

            deform_character_anim_bones(
            simulator, 
            curr_character_bone_anim_positions, 
            curr_character_bone_anim_rotations);
        }

        vel = inv_quat(curr_character_bone_anim_rotations(Bone_Entity)) 
                    * vec_to_vel(last_character_bone_anim_positions(Bone_Entity), curr_character_bone_anim_positions(Bone_Entity));

        avel = inv_quat(curr_character_bone_anim_rotations(Bone_Entity)) 
                    * quat_to_avel(last_character_bone_anim_rotations(Bone_Entity), curr_character_bone_anim_rotations(Bone_Entity));

        ltoe = inv_quat(curr_character_bone_anim_rotations(Bone_Entity)) 
                    * (curr_character_bone_anim_positions(Bone_LeftToe) - curr_character_bone_anim_positions(Bone_Entity));

        rtoe = inv_quat(curr_character_bone_anim_rotations(Bone_Entity)) 
                    * (curr_character_bone_anim_positions(Bone_RightToe) - curr_character_bone_anim_positions(Bone_Entity));

        if(t > 0)
        {
            lvel = vec_to_vel(features(t - 1, 0), ltoe);
            rvel = vec_to_vel(features(t - 1, 1), rtoe);

            // if(contact == 0 && dot(vel, rvel) < 0.f)
            // {
            //     contact = 1;
            //     contact_point = curr_character_bone_anim_positions(Bone_RightToe);
            // }
            // if(contact == 1 && dot(vel, lvel) < 0.f)
            // {
            //     contact = 0;
            //     contact_point = curr_character_bone_anim_positions(Bone_LeftToe);
            // }

            if(dot(vel, rvel) < dot(vel, lvel))
            {
                contact = 1;
                contact_point = curr_character_bone_anim_positions(Bone_RightToe);
            }
            else
            {
                contact = 0;
                contact_point = curr_character_bone_anim_positions(Bone_LeftToe);
            }
        }

        // if(contact == 0)
        // {
        //     IK_two_bones(
        //         // curr_character_bone_anim_positions(Bone_LeftUpLeg), 
        //         curr_character_bone_anim_positions(Bone_LeftLeg),
        //         curr_character_bone_anim_positions(Bone_LeftFoot), 
        //         curr_character_bone_anim_positions(Bone_LeftToe), 
        //         // curr_character_bone_anim_rotations(Bone_LeftUpLeg), 
        //         curr_character_bone_anim_rotations(Bone_LeftLeg), 
        //         curr_character_bone_anim_rotations(Bone_LeftFoot), 
        //         curr_character_bone_anim_rotations(Bone_LeftToe), 
        //         contact_point);
        // }
        // else
        // {
        //     IK_two_bones(
        //         // curr_character_bone_anim_positions(Bone_RightUpLeg),
        //         curr_character_bone_anim_positions(Bone_RightLeg),
        //         curr_character_bone_anim_positions(Bone_RightFoot), 
        //         curr_character_bone_anim_positions(Bone_RightToe), 
        //         // curr_character_bone_anim_rotations(Bone_RightUpLeg), 
        //         curr_character_bone_anim_rotations(Bone_RightLeg), 
        //         curr_character_bone_anim_rotations(Bone_RightFoot), 
        //         curr_character_bone_anim_rotations(Bone_RightToe), 
        //         contact_point);
        // }

        // batch_forward_kinematics_full(
        //     active_motion, 
        //     frame_num, 
        //     curr_character_bone_anim_positions, 
        //     curr_character_bone_anim_rotations);

        deform_character_anim_mesh(
            character, 
            curr_character_bone_anim_positions, 
            curr_character_bone_anim_rotations, 
            mesh1);

        // float delta = 0;
        for(int i = 2; i < 24; i++)
        {
            // renderer.models[i]->transform = mat4(eye3(), x(t, i - 2)) * mat4(R(t, i - 2), vec3());
            renderer.models[i]->transform = mat4(eye3(), simulator.bone_shapes(i - 2).pos) * mat4(quat_to_Rodrigues(simulator.bone_shapes(i - 2).rot), vec3());
            // delta = delta + length(x(t, i - 2) - simulator.bone_shapes(i - 2).pos);
        }
        // std::cout << t << ',' << delta << std::endl;

        // camera_controller.dir_gamepad_control(vec3(0.05f, 0.f, 0.f), 0.1f);
        camera_controller.dir_gamepad_control(gamepadstick_right, 0.1f);
        camera.update_orientation(vec3(0.f, 1.f, 0.f));
        // camera_controller.pos_pid_control(
        //     1.f, 
        //     0.01f, 
        //     0.1f, 
        //     500.f, 
        //     dt, 
        //     curr_character_bone_anim_positions(Bone_Hips) + camera.orientation * vec3(0.f, 0.f, 2.5f));
        camera_controller.pos_gamepad_control(vec3(curr_character_bone_anim_positions(Bone_Hips).x, 0.35f, curr_character_bone_anim_positions(Bone_Hips).z), camera.orientation, 2.5);

        light.light_position = curr_character_bone_anim_positions(Bone_Hips) + camera.orientation * vec3(0.f, 0.f, 5.f);
        light.light_radiance = 20.f + 10.f * sin(deg_to_rad(0.5f * t));

        renderer.models[0]->transform = mat4(eye3(), renderer.point_lights[0]->light_position) * mat4(0.01f);

        // draw(renderer.models[0], camera, &camera.fbo[0], &par, false);
        // draw(renderer.models[1], camera, &camera.fbo[0], &par, false);
        // for(int i = 2; i < 24; i++)
        // {
        //     draw(renderer.models[i], camera, &camera.fbo[0], &par, false);
        // }
        draw(renderer.models[0], camera, &camera.fbo[0], &par, true);
        draw(renderer.models[1], camera, &camera.fbo[0], &par, true);
        // for(int i = 2; i < 24; i++)
        // {
        //     draw(renderer.models[i], camera, &camera.fbo[0], &par, true);
        // }

        if(contact == 0)
        {
            draw_point(curr_character_bone_anim_positions(Bone_LeftUpLeg), camera, &camera.fbo[0], vec3(0.f, 255.f, 0.f));
            draw_point(curr_character_bone_anim_positions(Bone_LeftLeg), camera, &camera.fbo[0], vec3(0.f, 255.f, 0.f));
            draw_point(curr_character_bone_anim_positions(Bone_LeftFoot), camera, &camera.fbo[0], vec3(0.f, 255.f, 0.f));
            draw_point(curr_character_bone_anim_positions(Bone_LeftToe), camera, &camera.fbo[0], vec3(0.f, 255.f, 0.f));
            draw_point(contact_point, camera, &camera.fbo[0], vec3(255.f, 0.f, 0.f));
        }
        else
        {
            draw_point(curr_character_bone_anim_positions(Bone_RightUpLeg), camera, &camera.fbo[0], vec3(0.f, 255.f, 0.f));
            draw_point(curr_character_bone_anim_positions(Bone_RightLeg), camera, &camera.fbo[0], vec3(0.f, 255.f, 0.f));
            draw_point(curr_character_bone_anim_positions(Bone_RightFoot), camera, &camera.fbo[0], vec3(0.f, 255.f, 0.f));
            draw_point(curr_character_bone_anim_positions(Bone_RightToe), camera, &camera.fbo[0], vec3(0.f, 255.f, 0.f));
            draw_point(contact_point, camera, &camera.fbo[0], vec3(0.f, 0.f, 255.f));
        }

        draw_curve(curr_character_bone_anim_positions(Bone_Entity), 
                   curr_character_bone_anim_rotations(Bone_Entity) * vel, 
                   global_gamepad_vel, 
                   camera, 
                   &camera.fbo[0], 
                   vec3(0.f, 255.f, 0.f), 
                   N);

        for(int i = 0; i < motion.nbones(); i++)
        {
            draw_point(curr_character_bone_anim_positions(i), camera, &camera.fbo[0], vec3(0.f, 255.f, 0.f));
        }

        cv::imshow("MOOLAB", fbo_to_img(&camera.fbo[0]));

        // frames(t) = t;
        // contacts(t) = contact;

        // vels(t, 0) = vel;
        // avels(t, 0) = avel;

        // for(int j = 1; j < 6; j++)
        // {
        //     batch_forward_kinematics_full(
        //         active_motion, 
        //         frame_num + dN * j, 
        //         next_character_bone_anim_positions, 
        //         next_character_bone_anim_rotations);
            
        //     vels(t, j) = inv_quat(curr_character_bone_anim_rotations(Bone_Entity)) 
        //            * vec_to_vel(curr_character_bone_anim_positions(Bone_Entity), next_character_bone_anim_positions(Bone_Entity), dN * j * 0.0166667f);

        //     avels(t, j) = inv_quat(curr_character_bone_anim_rotations(Bone_Entity)) 
        //             * quat_to_avel(curr_character_bone_anim_rotations(Bone_Entity), next_character_bone_anim_rotations(Bone_Entity), dN * j * 0.0166667f);
        // }

        // features(t, 0) = ltoe;
        // features(t, 1) = rtoe;
        // features(t, 2) = lvel;
        // features(t, 3) = rvel;

        std::cout << " Time: " << t * dt << " Frame: " << frame_num 
                //   << " vel: " << db.vels(t, 0).x << ", " << db.vels(t, 0).y << ", " << db.vels(t, 0).z
                //   << " avel: " << db.avels(t, 0).x << ", " << db.avels(t, 0).y << ", " << db.avels(t, 0).z
        //           << " dir: " << camera.direction.x << ", " << camera.direction.y << ", " << camera.direction.z
        //           << " vel: " << vel.x << ", " << vel.y << ", " << vel.z 
        //           << " avel: " << avel.x << ", " << avel.y << ", " << avel.z 
        //           << " left: " << gamepadstick_left.x << ", " << gamepadstick_left.z 
        //           << " right: " << gamepadstick_right.x << ", " << gamepadstick_right.z
                  << std::endl;
                  
        frame_num++;
        t++;
        key = cv::waitKey(1);
        EndDrawing();
    }

    CloseWindow();

    // features(0, 0) = vec3(0.f);
    // features(0, 1) = vec3(0.f);
    // features(281, 0) = vec3(0.f);
    // features(281, 1) = vec3(0.f);
    // features(562, 0) = vec3(0.f);
    // features(562, 1) = vec3(0.f);
    // features(13153, 0) = vec3(0.f);
    // features(13153, 1) = vec3(0.f);
    // features(25744, 0) = vec3(0.f);
    // features(25744, 1) = vec3(0.f);
    // features(39622, 0) = vec3(0.f);
    // features(39622, 1) = vec3(0.f);

    // std::ofstream outfile;
    // outfile.open("../resources/motiondata.txt");

    // for(int t = 0; t < motion.nframes() - N; t++)
    // {
    //     outfile << features(t, 0).x << " " << features(t, 0).y << " " << features(t, 0).z << " " 
    //             << features(t, 1).x << " " << features(t, 1).y << " " << features(t, 1).z << " " 
    //             << features(t, 2).x << " " << features(t, 2).y << " " << features(t, 2).z << " "
    //             << features(t, 3).x << " " << features(t, 3).y << " " << features(t, 3).z << "\n";
    // }

    // outfile.close();

    // std::ifstream infile;
    // infile.open("../resources/output.txt");

    // for(int i = 0; i < motion.nframes() - N; i++)
    // {
    //     infile >> frames(i);
    //     infile >> velocities(i);
    //     // std::cout << frames(i) << ',' << velocities(i) << std::endl;
    // }

    // infile.close();

    // FILE* fft = fopen("../resources/features.bin", "rb");
    // assert(fft != NULL);
    // array2d_read(features, fft);
    // fclose(fft);

    // db.frames = frames;
    // db.contacts = contacts;
    // db.vels = vels;
    // db.avels = avels;
    // db.features = features;

    // std::cout << frames.size << ", " << velocities.size << ", " << features.rows << ", " << features.cols << std::endl;

    // Database_sort(db, features);

    // std::cout << db.frames.size << ", " << db.velocities.size << ", " << db.features.rows << ", " << db.features.cols << std::endl;

    // Database_revise(db);
    // Database_save(db, "../resources/database.bin");

    return 0;
}
