extern "C"
{
#include "raylib.h"
#include "raymath.h"
}
<<<<<<< HEAD

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

=======

#include "motion.hpp"
#include "character.hpp"
#include "curve.hpp"
#include "database.hpp"
#include "ik.hpp"
#include "motion.hpp"
#include "nnet.hpp"
#include "renderer.hpp"
#include "robot.hpp"

>>>>>>> 3795e8f (v1.04)
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
<<<<<<< HEAD

array1d<quat> transform_float_quat(const slice1d<float> res, vec3 axis, quat q)
{
    array1d<quat> arr(res.size);
    for(int i = 0; i < res.size; i++)
    {
        arr(i) = q;
        // arr(i) = quat(res(i), axis);
        // arr(i) = q * quat(res(i), axis);
        // std::cout << res(i) << std::endl;
    }
    return arr;
}
=======
>>>>>>> 3795e8f (v1.04)

int main(int argc, char** argv)
{
    InitWindow(720, 720, "raylib [background]");
    SetTargetFPS(60);
    //相机
<<<<<<< HEAD
    MooCamera camera(1, 0.1f, 0.1f, 0.05f, 50.f, vec3(0.f, 1.3f, -6.f), mat3());
    camera.set_view(vec3(0.f, 0.f, 1.f));
=======
    MooCamera camera(1, 0.1f, 0.1f, 0.05f, 50.f, vec3(0.f, 0.f, 3.f), mat3(), 720, 720);
    camera.set_view(vec3(0.f, 0.f, -1.f));
>>>>>>> 3795e8f (v1.04)
    camera.update_orientation(vec3(0.f, 1.f, 0.f));
    //控制器
    Camera_Controller camera_controller;
    bind_controller(camera_controller, &camera.position, &camera.direction);
<<<<<<< HEAD
    Character_Controller character_controller;
    character_controller.rotation0 = quat(1.f, 0.f, 0.f, 0.f);
    //角色，保存原始所有信息
    Character character;
    character_load(character, "../resources/character.bin");
    //数据库
    Database db;
    Database_load(db, "../resources/database.bin");
    //对角色数据的拷贝
    MooMesh mesh0 = make_character_rest_mesh(character);
    MooMesh mesh1 = cube();
    //材质
    MooMaterial material0(PBR);
    material0.roughness = 0.1f;
    material0.tex = img_to_tex("../resources/skybox.png");
    material0.BRDFLut = img_to_tex("../resources/BRDFLut.png");
    material0.EavgLut = img_to_tex("../resources/EavgLut.png");
    MooMaterial material1(TEXTURE);
    material1.roughness = 0.1f;
    material1.tex = img_to_tex("../resources/UI/0.png");
    material1.BRDFLut = img_to_tex("../resources/BRDFLut.png");
    material1.EavgLut = img_to_tex("../resources/EavgLut.png");
    //model存储mesh的拷贝，但由于mesh本身存储的是地址，实际上仍为mesh地址上的数据
    MooModel model1 = mesh_material_to_model(mesh0, material0);
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
    // active_motion = motion;
    search_motion = motion_sub_sequence(motion, 0, 0 + 1.5 * N);
    // search_motion = motion_sub_sequence(motion, 0, 281);
    // search_motion = motion_sub_sequence(motion, 281, 562);
    // search_motion = motion_sub_sequence(motion, 562, 13153);
    // search_motion= motion_sub_sequence(motion, 13153, 25744);
    // search_motion = motion_sub_sequence(motion, 25744, 39622);
    // search_motion = motion_sub_sequence(motion, 39622, 53500);
    // search_motion = build_loop_motion(active_motion, 0.1f, 0.1f);
    //模拟器
    Simulator simulator;
    simulator.bind_simulator(motion);
    //点光源
    Point_Light light(vec3(0.f, 1.f, 1.5f), vec3(0.f, 0.f, -1.f), vec3(20.f));
    MooModel model0 = get_light_model(light);
    MooModel model2 = mesh_material_to_model(mesh1, material1);
    model2.transform = mat4(eye3(), vec3(1.5f, 0.f, 1.f)) * mat4(0.25f, 0.25f, 0.25f);
=======
    //点光源
    Point_Light light(vec3(0.f, 1.f, 1.5f), vec3(0.f, 0.f, -1.f), vec3(20.f));
    //对角色数据的拷贝
    MooMesh mesh1 = capsule(0.5f * 0.0738f, 0.0738f, 18, 18,  0.5f * 0.0738f);
    MooMesh mesh2 = capsule(0.5f * 0.1f, 0.435f, 18, 18,  0.5f * 0.1f);
    MooMesh mesh3 = capsule(0.5f * 0.1f, 0.4237f, 18, 18,  0.5f * 0.1f);
    MooMesh mesh4 = capsule(0.5f * 0.1f, 0.1730f, 18, 18,  0.5f * 0.1f);
    MooMesh mesh5 = capsule(0.5f * 0.05f, 0.05f, 18, 18,  0.5f * 0.05f);
    MooMesh mesh6 = capsule(0.5f * 0.12f, 0.1259f, 18, 18,  0.5f * 0.12f);
    MooMesh mesh7 = capsule(0.5f * 0.12f, 0.1234f, 18, 18,  0.5f * 0.12f);
    MooMesh mesh8 = capsule(0.5f * 0.12f, 0.2583f, 18, 18,  0.5f * 0.12f);
    MooMesh mesh9 = capsule(0.5f * 0.10f, 0.1177f, 18, 18,  0.5f * 0.10f);
    MooMesh mesh10 = capsule(0.5f * 0.15f, 0.25f, 18, 18,  0.5f * 0.15f);
    MooMesh mesh11 = capsule(0.5f * 0.1128f, 0.2583f, 18, 18,  0.5f * 0.1128f);
    MooMesh mesh12 = capsule(0.5f * 0.1f, 0.33f, 18, 18,  0.5f * 0.1f);
    MooMesh mesh13 = capsule(0.5f * 0.1f, 0.252f, 18, 18,  0.5f * 0.1f);
    MooMesh mesh14 = capsule(0.5f * 0.1f, 0.2f, 18, 18,  0.5f * 0.1f);
    //材质
    MooMaterial material1(TEXTURE);
    material1.roughness = 0.1f;
    material1.tex = img_to_tex("../resources/skybox.png");
    material1.BRDFLut = img_to_tex("../resources/BRDFLut.png");
    material1.EavgLut = img_to_tex("../resources/EavgLut.png");
    //model存储mesh的拷贝，但由于mesh本身存储的是地址，实际上仍为mesh地址上的数据
    MooModel model0 = get_light_model(light);
    MooModel model1 = mesh_material_to_model(mesh1, material1);
    MooModel model2 = mesh_material_to_model(mesh2, material1);
    MooModel model3 = mesh_material_to_model(mesh3, material1);
    MooModel model4 = mesh_material_to_model(mesh4, material1);
    MooModel model5 = mesh_material_to_model(mesh5, material1);
    MooModel model6 = mesh_material_to_model(mesh2, material1);
    MooModel model7 = mesh_material_to_model(mesh3, material1);
    MooModel model8 = mesh_material_to_model(mesh4, material1);
    MooModel model9 = mesh_material_to_model(mesh5, material1);
    MooModel model10 = mesh_material_to_model(mesh6, material1);
    MooModel model11 = mesh_material_to_model(mesh7, material1);
    MooModel model12 = mesh_material_to_model(mesh8, material1);
    MooModel model13 = mesh_material_to_model(mesh9, material1);
    MooModel model14 = mesh_material_to_model(mesh10, material1);
    MooModel model15 = mesh_material_to_model(mesh11, material1);
    MooModel model16 = mesh_material_to_model(mesh12, material1);
    MooModel model17 = mesh_material_to_model(mesh13, material1);
    MooModel model18 = mesh_material_to_model(mesh14, material1);
    MooModel model19 = mesh_material_to_model(mesh11, material1);
    MooModel model20 = mesh_material_to_model(mesh12, material1);
    MooModel model21 = mesh_material_to_model(mesh13, material1);
    MooModel model22 = mesh_material_to_model(mesh14, material1);
>>>>>>> 3795e8f (v1.04)
    //更新信息
    updated_paramters par;
    par.light_dir = &light.light_direction;
    par.light_pos = &light.light_position;
    par.light_radiance = &light.light_radiance;
    par.shadow_map = light.shadow_map;
<<<<<<< HEAD
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
=======
>>>>>>> 3795e8f (v1.04)
    //渲染器
    MooRenderer renderer;

    renderer.models.push_back(&model0);
    renderer.models.push_back(&model1);
    renderer.models.push_back(&model2);
<<<<<<< HEAD
    renderer.point_lights.push_back(&light);

    FBO output0(720, 720);
=======
    renderer.models.push_back(&model3);
    renderer.models.push_back(&model4);
    renderer.models.push_back(&model5);
    renderer.models.push_back(&model6);
    renderer.models.push_back(&model7);
    renderer.models.push_back(&model8);
    renderer.models.push_back(&model9);
    renderer.models.push_back(&model10);
    renderer.models.push_back(&model11);
    renderer.models.push_back(&model12);
    renderer.models.push_back(&model13);
    renderer.models.push_back(&model14);
    renderer.models.push_back(&model15);
    renderer.models.push_back(&model16);
    renderer.models.push_back(&model17);
    renderer.models.push_back(&model18);
    renderer.models.push_back(&model19);
    renderer.models.push_back(&model20);
    renderer.models.push_back(&model21);
    renderer.models.push_back(&model22);
    renderer.point_lights.push_back(&light);
>>>>>>> 3795e8f (v1.04)

    int key = 0;
    int frame_num = 0;
    int t = 0;

<<<<<<< HEAD
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

    int weight_vel = 10;
    int weight_avel = 6;
    int weight_toe = 3;

    batch_forward_kinematics_full(
        motion, 
        0, 
        curr_character_bone_anim_positions, 
        curr_character_bone_anim_rotations);
    contact = 0;
    contact_point = curr_character_bone_anim_positions(Bone_LeftToe);

    array1d<float> alpha(0.5 * N);
    for(int t = 0; t < 0.5 * N; t++)
    {
        alpha(t) = float(t) / (0.5 * N);
    }

    float dt = 0.0166667f;

    // std::cout << db.frames.size << ',' << db.velocities.size << ',' << db.features.rows << ',' << db.features.cols << std::endl;

    array1d<float> sh_points(5);
    array1d<float> el_points(5);
    array1d<float> wr_points(5);
    array1d<float> td(4);

    sh_points(0) = 105.f;
    sh_points(1) = 60.f;
    sh_points(2) = 15.f;
    sh_points(3) = 60.f;
    sh_points(4) = 105.f;

    el_points(0) = -80.f;
    el_points(1) = -100.f;
    el_points(2) = 0.f;
    el_points(3) = -70.f;
    el_points(4) = -80.f;

    wr_points(0) = 0.f;
    wr_points(1) = -60.f;
    wr_points(2) = 90.f;
    wr_points(3) = -10.f;
    wr_points(4) = 0.f;

    td(0) = 0.4f;
    td(1) = 0.1f;
    td(2) = 0.8f;
    td(3) = 1.2f;

    array1d<float> sh_curve_points = LFPB(sh_points, td, 10000, 0.0166667f);
    array1d<float> el_curve_points = LFPB(el_points, td, 10000, 0.0166667f);
    array1d<float> wr_curve_points = LFPB(wr_points, td, 10000, 0.0166667f);

    array1d<quat> sh_curve_quaternions = transform_float_quat(sh_curve_points, vec3(0.f, 0.f, 1.f), quat(90.f, vec3(1.f, 0.f, 0.f)) * quat(-90.f, vec3(0.f, 0.f, 1.f)) * quat(90.f, vec3(1.f, 0.f, 0.f)));
    array1d<quat> el_curve_quaternions = transform_float_quat(el_curve_points, vec3(0.f, 1.f, 0.f), quat(1.f, 0.f, 0.f, 0.f));
    array1d<quat> wr_curve_quaternions = transform_float_quat(wr_curve_points, vec3(0.f, 1.f, 0.f), quat(-90.f, vec3(1.f, 0.f, 0.f)));

    // active_motion = motion;

    cv::imshow("MOOLAB", fbo_to_img(&output0));
    cv::createTrackbar("vel", "MOOLAB", &weight_vel, 50);
    cv::createTrackbar("avel", "MOOLAB", &weight_avel, 50);
    cv::createTrackbar("toe", "MOOLAB", &weight_toe, 50);

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

        // if(frame_num == N)
        // {
            // vel = (2.5f + 2.5f * sin(deg_to_rad(0.1f * t))) * vec3(0.f, 0.f, 1.f);

            // avel = vec3(0.f, 0.f, 0.f);

            // std::cout << 2.5f + 2.5f * sin(deg_to_rad(0.1f * t)) << ',' << vel.x << ',' << vel.y << ',' << vel.z << std::endl;
        // }
        
        if(frame_num == N)
        {
            frame_num = 0;

            // character_controller.vel5 = vel;
            // character_controller.avel5 = avel;
            character_controller.update(local_gamepad_vel, local_gamepad_avel, curr_character_bone_anim_rotations(Bone_Entity));
            int best_frame = Database_search(db, character_controller, contact, ltoe, rtoe, N, dN, weight_vel, weight_avel, weight_toe);
            search_motion = motion_sub_sequence(
                                    motion, 
                                    best_frame, 
                                    best_frame + 1.5 * N);
            active_motion = motion_sub_sequence(
                                concatenate_two_motions(
                                    active_motion, 
                                    search_motion, 
                                    N, 
                                    0.5 * N), 
                                N,
                                N + 1.5 * N);
            // active_motion = search_motion;
            // array2d_set_by_cols(active_motion.bone_local_rotations, sh_curve_quaternions, Bone_RightArm);
            // array2d_set_by_cols(active_motion.bone_local_rotations, el_curve_quaternions, Bone_RightForeArm);
            // array2d_set_by_cols(active_motion.bone_local_rotations, wr_curve_quaternions, Bone_RightHand);
        }

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

        output0.set(vec3(100.f, 100.f, 200.f), 100.f);

        last_character_bone_anim_positions = curr_character_bone_anim_positions;
        last_character_bone_anim_rotations = curr_character_bone_anim_rotations;

        batch_forward_kinematics_full(
            active_motion, 
            frame_num, 
            curr_character_bone_anim_positions, 
            curr_character_bone_anim_rotations);

        // batch_forward_kinematics_part(
        //     active_motion, 
        //     frame_num, 
        //     curr_character_bone_anim_positions, 
        //     curr_character_bone_anim_rotations, 
        //     vec3(0.f), 
        //     quat(0.f, vec3(0.f, 0.f, 1.f)));

        // batch_forward_kinematics_root(
        //     active_motion, 
        //     frame_num, 
        //     curr_character_bone_anim_positions, 
        //     curr_character_bone_anim_rotations);

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

        if(contact == 0)
        {
            IK_two_bones(
                // curr_character_bone_anim_positions(Bone_LeftUpLeg), 
                curr_character_bone_anim_positions(Bone_LeftLeg),
                curr_character_bone_anim_positions(Bone_LeftFoot), 
                curr_character_bone_anim_positions(Bone_LeftToe), 
                // curr_character_bone_anim_rotations(Bone_LeftUpLeg), 
                curr_character_bone_anim_rotations(Bone_LeftLeg), 
                curr_character_bone_anim_rotations(Bone_LeftFoot), 
                curr_character_bone_anim_rotations(Bone_LeftToe), 
                contact_point);
        }
        else
        {
            IK_two_bones(
                // curr_character_bone_anim_positions(Bone_RightUpLeg),
                curr_character_bone_anim_positions(Bone_RightLeg),
                curr_character_bone_anim_positions(Bone_RightFoot), 
                curr_character_bone_anim_positions(Bone_RightToe), 
                // curr_character_bone_anim_rotations(Bone_RightUpLeg), 
                curr_character_bone_anim_rotations(Bone_RightLeg), 
                curr_character_bone_anim_rotations(Bone_RightFoot), 
                curr_character_bone_anim_rotations(Bone_RightToe), 
                contact_point);
        }

        // simulator.simulate_gravity();
        // simulator.simulate();
        // simulator.batch_forward_kinematics_full();
        // curr_character_bone_anim_positions = simulator.bone_anim_positions; 
        // curr_character_bone_anim_rotations = simulator.bone_anim_rotations; 

        deform_character_anim_mesh(
            character, 
            curr_character_bone_anim_positions, 
            curr_character_bone_anim_rotations, 
            mesh0);

        // light.light_position = curr_character_bone_anim_positions(Bone_Entity) + 
        //                        curr_character_bone_anim_rotations(Bone_Entity) * vec3(0.f, 1.5f, -1.5f);
        // light.light_radiance = 20.f + 10.f * sin(deg_to_rad(0.5f * t));

        // camera_controller.pos_pid_control(
        //     1.f, 
        //     0.01f, 
        //     0.1f, 
        //     500.f, 
        //     dt, 
        //     curr_character_bone_anim_positions(Bone_Entity) + curr_character_bone_anim_rotations(Bone_Entity) * vec3(0.f, 1.2f, -2.f));
        // camera_controller.dir_pid_control(
        //     0.7f, 
        //     0.01f, 
        //     0.1f, 
        //     500.f, 
        //     dt, 
        //     curr_character_bone_anim_rotations(Bone_Entity) * vec3(0.f, 0.f, 1.f));
        // camera.update_orientation(vec3(0.f, 1.f, 0.f));

        // camera.position = curr_character_bone_anim_positions(Bone_Entity) + curr_character_bone_anim_rotations(Bone_Entity) * vec3(0.f, 1.2f, -2.f);
        // camera.direction = curr_character_bone_anim_rotations(Bone_Entity) * vec3(0.f, 0.f, 1.f);
        // camera.update_orientation(vec3(0.f, 1.f, 0.f));

        camera_controller.dir_gamepad_control(gamepadstick_right, 0.1f);
        camera.update_orientation(vec3(0.f, 1.f, 0.f));
        // camera.position = curr_character_bone_anim_positions(Bone_Hips) + camera.orientation * vec3(0.f, 0.f, 2.f);
        camera_controller.pos_pid_control(
            1.f, 
            0.01f, 
            0.1f, 
            500.f, 
            dt, 
            curr_character_bone_anim_positions(Bone_Hips) + camera.orientation * vec3(0.f, 0.f, 2.5f));

        light.light_position = curr_character_bone_anim_positions(Bone_Hips) + camera.orientation * vec3(0.f, 0.f, 5.f);

        renderer.models[0]->transform = mat4(eye3(), renderer.point_lights[0]->light_position) * mat4(0.01f);
        // renderer.models[0]->transform = mat4(eye3(), curr_character_bone_anim_positions(Bone_Entity)) * mat4(0.1f);

        draw(renderer.models[0], camera, &output0, &par, true);
        draw(renderer.models[1], camera, &output0, &par, true);
        // draw(renderer.models[2], camera, &output0, &par, true);

        if(contact == 0)
        {
            draw_point(curr_character_bone_anim_positions(Bone_LeftUpLeg), camera, &output0, vec3(0.f, 255.f, 0.f));
            draw_point(curr_character_bone_anim_positions(Bone_LeftLeg), camera, &output0, vec3(0.f, 255.f, 0.f));
            draw_point(curr_character_bone_anim_positions(Bone_LeftFoot), camera, &output0, vec3(0.f, 255.f, 0.f));
            draw_point(curr_character_bone_anim_positions(Bone_LeftToe), camera, &output0, vec3(0.f, 255.f, 0.f));
            draw_point(contact_point, camera, &output0, vec3(255.f, 0.f, 0.f));
        }
        else
        {
            draw_point(curr_character_bone_anim_positions(Bone_RightUpLeg), camera, &output0, vec3(0.f, 255.f, 0.f));
            draw_point(curr_character_bone_anim_positions(Bone_RightLeg), camera, &output0, vec3(0.f, 255.f, 0.f));
            draw_point(curr_character_bone_anim_positions(Bone_RightFoot), camera, &output0, vec3(0.f, 255.f, 0.f));
            draw_point(curr_character_bone_anim_positions(Bone_RightToe), camera, &output0, vec3(0.f, 255.f, 0.f));
            draw_point(contact_point, camera, &output0, vec3(0.f, 0.f, 255.f));
        }

        draw_curve(curr_character_bone_anim_positions(Bone_Entity), 
                   curr_character_bone_anim_rotations(Bone_Entity) * vel, 
                   global_gamepad_vel, 
                   camera, 
                   &output0, 
                   vec3(0.f, 255.f, 0.f), 
                   N);

        // UI_draw_char('4', vec2(360, 360), ui, &output0, vec3(0.f, 255.f, 0.f));

        cv::imshow("MOOLAB", fbo_to_img(&output0));

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

        features(t, 0) = ltoe;
        features(t, 1) = rtoe;
        features(t, 2) = lvel;
        features(t, 3) = rvel;

        std::cout << " Time: " << t << " Frame: " << frame_num 
                  << " vel: " << db.vels(t, 0).x << ", " << db.vels(t, 0).y << ", " << db.vels(t, 0).z
                  << " avel: " << db.avels(t, 0).x << ", " << db.avels(t, 0).y << ", " << db.avels(t, 0).z
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
=======
    float dt = 0.0166667f;

    int times = 3000;

    array2d<vec3> x(times, 22);
    array2d<mat3> R(times, 22);

    std::ifstream infile;

    infile.open("../resources/x.txt");
    for(int i = 0; i < times; i++)
    {
        for(int j = 0; j < 22; j++)
        {
            infile >> x(i, j).x;
        }
        for(int j = 0; j < 22; j++)
        {
            infile >> x(i, j).y;
        }
        for(int j = 0; j < 22; j++)
        {
            infile >> x(i, j).z;
        }
    }
    infile.close();

    infile.open("../resources/R.txt");
    for(int i = 0; i < times; i++)
    {
        for(int j = 0; j < 22; j++)
        {
            infile >> R(i, j).X.x;
            infile >> R(i, j).Y.x;
            infile >> R(i, j).Z.x;
        }
        for(int j = 0; j < 22; j++)
        {
            infile >> R(i, j).X.y;
            infile >> R(i, j).Y.y;
            infile >> R(i, j).Z.y;
        }
        for(int j = 0; j < 22; j++)
        {
            infile >> R(i, j).X.z;
            infile >> R(i, j).Y.z;
            infile >> R(i, j).Z.z;
        }
    }
    infile.close();

    cv::imshow("MOOLAB", fbo_to_img(&camera.fbo[0]));

    while (key != 27 && t < times)
    {
        vec3 gamepadstick_left = gamepad_get_stick(GAMEPAD_STICK_LEFT);
        vec3 gamepadstick_right = gamepad_get_stick(GAMEPAD_STICK_RIGHT);

        BeginDrawing();
        ClearBackground(RAYWHITE);

        camera.fbo[0].set(vec3(100.f, 100.f, 200.f), 100.f);

        // camera.position = camera.position + 5.f * (quat(rad_to_deg(camera_controller.ang.x), vec3(0.f, 1.f, 0.f)) * (-gamepadstick_left)) * dt;
        // camera_controller.dir_gamepad_control(gamepadstick_right, 0.1f);
        // camera.update_orientation(vec3(0.f, 1.f, 0.f));
        // camera.position = vec3(0.f, x(t, 0).y, 3.f);
        camera_controller.pos_pid_control(1.f, 0.05f, 0.1f, 500.f, dt, vec3(x(t, 0).x, x(t, 0).y, x(t, 0).z + 3.f));

        renderer.models[0]->transform = mat4(eye3(), renderer.point_lights[0]->light_position) * mat4(0.01f);
        for(int i = 1; i < 23; i++)
        {
            if((i >= 2 && i <= 5) || (i >= 15 && i <= 18))
            {
                renderer.models[i]->transform = mat4(Rodrigues(0, vec3(0.f, 1.f, 0.f)), vec3()) *  mat4(eye3(), x(t, i - 1)) * mat4(R(t, i - 1), vec3());
            }
            else if((i >= 6 && i <= 9) || (i >= 19 && i <= 22))
            {
                renderer.models[i]->transform = mat4(Rodrigues(0, vec3(0.f, 1.f, 0.f)), vec3()) *  mat4(eye3(), x(t, i - 1)) * mat4(R(t, i - 1), vec3());
            }
            else
            {
                renderer.models[i]->transform = mat4(Rodrigues(0, vec3(0.f, 1.f, 0.f)), vec3()) *  mat4(eye3(), x(t, i - 1)) * mat4(R(t, i - 1), vec3());
            }
        }
        
        draw(renderer.models[0], camera, &camera.fbo[0], &par, false);
        for(int i = 1; i < 23; i++)
        {
            draw(renderer.models[i], camera, &camera.fbo[0], &par, false);
        }
        draw(renderer.models[0], camera, &camera.fbo[0], &par, true);
        for(int i = 1; i < 23; i++)
        {
            draw(renderer.models[i], camera, &camera.fbo[0], &par, true);
        }

        cv::imshow("MOOLAB", fbo_to_img(&camera.fbo[0]));

        std::cout << " Time: " << t << " Frame: " << frame_num 
                  << std::endl;
        frame_num++;
        t++;
        key = cv::waitKey(1);
        EndDrawing();
    }

    CloseWindow();
>>>>>>> 3795e8f (v1.04)

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
