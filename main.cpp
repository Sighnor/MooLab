extern "C"
{
#include "raylib.h"
#include "raymath.h"
}

#include "motion.hpp"
#include "character.hpp"
#include "controller.hpp"
#include "curve.hpp"
#include "database.hpp"
#include "motion.hpp"
#include "nnet.hpp"
#include "renderer.hpp"
#include "robot.hpp"

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
    float gamepadmag = sqrtf(gamepadx*gamepadx + gamepady*gamepady);
    
    if (gamepadmag > deadzone)
    {
        float gamepaddirx = gamepadx / gamepadmag;
        float gamepaddiry = gamepady / gamepadmag;
        float gamepadclippedmag = gamepadmag > 1.0f ? 1.0f : gamepadmag*gamepadmag;
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

array1d<quat> transform_float_quat(const array1d<float> res, vec3 axis, quat q)
{
    array1d<quat> arr(res.size);
    for(int i = 0; i < res.size; i++)
    {
        // arr(i) = q;
        // arr(i) = quat(res(i), axis);
        arr(i) = quat(res(i), axis) * q;
        // std::cout << res(i) << std::endl;
    }
    return arr;
}

int main(int argc, char** argv)
{
    InitWindow(720, 720, "raylib [background]");
    SetTargetFPS(60);
    //相机
    MooCamera camera(1, 0.1f, 0.1f, 0.05f, 50.f, vec3(0.f, 1.3f, -6.f), mat3());
    camera.set_view(vec3(0.f, 0.f, 1.f));
    camera.update_orientation(vec3(0.f, 1.f, 0.f));
    //控制器
    Controller controller;
    bind_controller(controller, &camera.position, &camera.direction);
    //角色，保存原始所有信息
    Character character;
    character_load(character, "../resources/character.bin");
    //数据库
    Database db;
    Database_load(db, "../resources/database.bin");
    //对角色数据的拷贝
    MooMesh mesh0 = make_character_rest_mesh(character);
    //材质
    MooMaterial material0(PBR);
    material0.roughness = 0.1f;
    material0.tex = img_to_tex("../resources/skybox.png");
    material0.BRDFLut = img_to_tex("../resources/BRDFLut.png");
    material0.EavgLut = img_to_tex("../resources/EavgLut.png");
    //model存储mesh的拷贝，但由于mesh本身存储的是地址，实际上仍为mesh地址上的数据
    MooModel model1 = mesh_material_to_model(mesh0, material0);
    //动作库
    int N = 100;
    BVH_Motion motion, motion1, motion2, motion3, motion4, motion5, motion6, active_motion, search_motion;
    Motion_load(motion, "../resources/long_motion.bin");
    motion1 = motion_sub_sequence(motion, 0, 250);
    motion2 = motion_sub_sequence(motion, 700, 1000);
    motion3 = translation_and_rotation(motion1, 0, vec3(0.f, 0.f, 0.f), vec3(0.f, 0.f, -1.f));
    motion4 = translation_and_rotation(motion2, 0, vec3(0.f, 0.f, 0.f), vec3(-1.f, 0.f, 0.f));
    motion5 = concatenate_two_motions(motion3, motion4, 150, 100);
    motion6 = motion;
    active_motion = translation_and_rotation(motion_sub_sequence(motion, 900, 900 + 1.5 * N), 0, vec3(-1.f, 0.f, 0.f), vec3(1.f, 0.f, 0.f));
    search_motion = translation_and_rotation(motion_sub_sequence(motion, 900, 900 + 1.5 * N), 0, vec3(-1.f, 0.f, 0.f), vec3(1.f, 0.f, 0.f));
    // motion6 = motion_sub_sequence(motion, 0, 281);
    // motion6 = motion_sub_sequence(motion, 281, 562);
    // motion6 = motion_sub_sequence(motion, 562, 13153);
    // motion6 = motion_sub_sequence(motion, 13153, 25744);
    // motion6 = motion_sub_sequence(motion, 25744, 39622);
    // motion6 = motion_sub_sequence(motion, 39622, 53500);
    //点光源
    Point_Light light(vec3(0.f, 1.f, 1.5f), vec3(0.f, 0.f, -1.f), vec3(20.f));
    MooModel model0 = get_light_model(light);
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
    renderer.point_lights.push_back(&light);

    FBO output0(720, 720);

    int key = 0;
    int frame_num = 0;
    int t = 0;

    array1d<vec3> curr_character_bone_anim_positions(motion.nbones());
    array1d<quat> curr_character_bone_anim_rotations(motion.nbones());

    array1d<vec3> last_character_bone_anim_positions(motion.nbones());
    array1d<quat> last_character_bone_anim_rotations(motion.nbones());
    // vel avel ltoe rtoe
    array1d<int> frames(motion.nframes());
    array1d<float> velocities(motion.nframes());
    array2d<vec3> features(motion.nframes(), 4);
    vec3 vel;
    vec3 avel;
    vec3 ltoe;
    vec3 rtoe;
    float contact = 0;
    vec3 contact_point;

    batch_forward_kinematics_full(
        motion, 
        0, 
        curr_character_bone_anim_positions, 
        curr_character_bone_anim_rotations);
    contact = 0;
    contact_point = curr_character_bone_anim_positions(Bone_LeftFoot);

    array1d<float> alpha(0.5 * N);
    for(int t = 0; t < 0.5 * N; t++)
    {
        alpha(t) = float(t) / (0.5 * N);
    }

    float dt = 1.f / 60.f;

    // std::cout << db.frames.size << ',' << db.velocities.size << ',' << db.features.rows << ',' << db.features.cols << std::endl;

    // array1d<float> sh_points(5);
    // array1d<float> el_points(5);
    // array1d<float> wr_points(5);
    // array1d<float> td(4);

    // sh_points(0) = 105.f;
    // sh_points(1) = 60.f;
    // sh_points(2) = 15.f;
    // sh_points(3) = 60.f;
    // sh_points(4) = 105.f;

    // el_points(0) = -80.f;
    // el_points(1) = -100.f;
    // el_points(2) = 0.f;
    // el_points(3) = -70.f;
    // el_points(4) = -80.f;

    // wr_points(0) = 0.f;
    // wr_points(1) = -60.f;
    // wr_points(2) = 90.f;
    // wr_points(3) = -10.f;
    // wr_points(4) = 0.f;

    // td(0) = 0.4f;
    // td(1) = 0.1f;
    // td(2) = 0.8f;
    // td(3) = 1.2f;

    // array1d<float> sh_curve_points = LFPB(sh_points, td, 10000, 0.0166667f);
    // array1d<float> el_curve_points = LFPB(el_points, td, 10000, 0.0166667f);
    // array1d<float> wr_curve_points = LFPB(wr_points, td, 10000, 0.0166667f);

    // array1d<quat> sh_curve_quaternions = transform_float_quat(sh_curve_points, vec3(1.f, 0.f, 0.f), quat(0.f, vec3(0.f, 1.f, 0.f)));
    // array1d<quat> el_curve_quaternions = transform_float_quat(el_curve_points, vec3(0.f, 0.f, 1.f), quat(1.f, 0.f, 0.f, 0.f));
    // array1d<quat> wr_curve_quaternions = transform_float_quat(wr_curve_points, vec3(0.f, 1.f, 0.f), quat(90.f, vec3(1.f, 0.f, 0.f)));

    while (key != 27 && t < motion.nframes())
    {
        vec3 gamepadstick_left = gamepad_get_stick(GAMEPAD_STICK_LEFT);
        vec3 gamepadstick_right = gamepad_get_stick(GAMEPAD_STICK_RIGHT);

        vel = 2.f * (inv_quat(curr_character_bone_anim_rotations(Bone_Entity)) 
                   * quat(rad_to_deg(controller.ang.x), vec3(0.f, 1.f, 0.f)) * (-gamepadstick_left));

        avel = vec3(0.f, 0.f, 0.f);

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
            int best_frame = Database_search(db, vel, avel, ltoe, rtoe);
            search_motion = translation_and_rotation(
                                motion_sub_sequence(
                                    motion, 
                                    best_frame, 
                                    best_frame + 1.5 * N), 
                                0, 
                                curr_character_bone_anim_positions(Bone_Entity) + vel * dt, 
                                quat(length(avel) * dt, normalize(avel)) * curr_character_bone_anim_rotations(Bone_Entity) * vec3(0.f, 0.f, 1.f));
            active_motion = motion_sub_sequence(
                                concatenate_two_motions(
                                    active_motion, 
                                    search_motion, 
                                    N, 
                                    0.1 * N), 
                                N,
                                N + 1.5 * N);
        }

        output0.set(vec3(100.f, 100.f, 200.f), 100.f);

        last_character_bone_anim_positions = curr_character_bone_anim_positions;
        last_character_bone_anim_rotations = curr_character_bone_anim_rotations;

        // batch_forward_kinematics_full(
        //     active_motion, 
        //     frame_num, 
        //     curr_character_bone_anim_positions, 
        //     curr_character_bone_anim_rotations);

        batch_forward_kinematics_part(
            active_motion, 
            frame_num, 
            curr_character_bone_anim_positions, 
            curr_character_bone_anim_rotations, 
            vec3(0.f), 
            quat(0.f, vec3(0.f, 0.f, 1.f)));

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

        // batch_forward_kinematics_root(
        //     active_motion, 
        //     frame_num, 
        //     curr_character_bone_anim_positions, 
        //     curr_character_bone_anim_rotations);

        deform_character_anim_mesh(
            character, 
            curr_character_bone_anim_positions, 
            curr_character_bone_anim_rotations, 
            mesh0);

        // vel = inv_quat(curr_character_bone_anim_rotations(Bone_Entity)) 
        //            * vec_to_vel(last_character_bone_anim_positions(Bone_Entity), curr_character_bone_anim_positions(Bone_Entity));

        // avel = inv_quat(curr_character_bone_anim_rotations(Bone_Entity)) 
        //             * quat_to_avel(last_character_bone_anim_rotations(Bone_Entity), curr_character_bone_anim_rotations(Bone_Entity));

        ltoe = inv_quat(curr_character_bone_anim_rotations(Bone_Entity)) 
                    * (curr_character_bone_anim_positions(Bone_LeftToe) - curr_character_bone_anim_positions(Bone_Entity));

        rtoe = inv_quat(curr_character_bone_anim_rotations(Bone_Entity)) 
                    * (curr_character_bone_anim_positions(Bone_RightToe) - curr_character_bone_anim_positions(Bone_Entity));

        // if(t > 0)
        // {
        //     vec3 lvel = vec_to_vel(features(t - 1, 2), ltoe);
        //     vec3 rvel = vec_to_vel(features(t - 1, 3), rtoe);
        //     if(contact == 0 && dot(vel, rvel) < 0.f)
        //     {
        //         contact = 1;
        //         contact_point = rtoe;
        //     }
        //     if(contact == 1 && dot(vel, lvel) < 0.f)
        //     {
        //         contact = 0;
        //         contact_point = ltoe;
        //     }
        // }

        if(t > 1 && ltoe.y > 0.05f || rtoe.y > 0.05f)
        {
            if(contact == 1 && ltoe.y < 0.05f)
            {
                contact = 0;
                contact_point = ltoe;
            }
            else if(contact == 0 && rtoe.y < 0.05f)
            {
                contact = 1;
                contact_point = rtoe;
            }
        }

        // light.light_position = curr_character_bone_anim_positions(Bone_Entity) + 
        //                        curr_character_bone_anim_rotations(Bone_Entity) * vec3(0.f, 1.5f, -1.5f);
        // light.light_radiance = 20.f + 10.f * sin(deg_to_rad(0.5f * frame_num));

        // controller.pos_pid_control(
        //     1.f, 
        //     0.01f, 
        //     0.1f, 
        //     500.f, 
        //     dt, 
        //     curr_character_bone_anim_positions(Bone_Entity) + curr_character_bone_anim_rotations(Bone_Entity) * vec3(0.f, 1.2f, -2.f));
        // controller.dir_pid_control(
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
        controller.dir_gamepad_control(gamepadstick_right, 0.1f);
        camera.update_orientation(vec3(0.f, 1.f, 0.f));
        camera.position = curr_character_bone_anim_positions(Bone_Hips) + camera.orientation * vec3(0.f, 0.f, 2.f);

        light.light_position = curr_character_bone_anim_positions(Bone_Hips) + camera.orientation * vec3(0.f, 0.f, 3.f);

        renderer.models[0]->transform = mat4(eye3(), renderer.point_lights[0]->light_position) * mat4(0.01f);
        // renderer.models[0]->transform = mat4(eye3(), curr_character_bone_anim_positions(Bone_Entity)) * mat4(0.1f);

        draw(renderer.models[0], camera, &output0, &par, true);
        draw(renderer.models[1], camera, &output0, &par, true);

        if(contact == 0)
        {
            draw_point(curr_character_bone_anim_positions(Bone_LeftUpLeg), camera, &output0, vec3(0.f, 255.f, 0.f));
            draw_point(curr_character_bone_anim_positions(Bone_LeftFoot), camera, &output0, vec3(0.f, 255.f, 0.f));
            draw_point(curr_character_bone_anim_positions(Bone_LeftToe), camera, &output0, vec3(0.f, 255.f, 0.f));
            // draw_point(contact_point, camera, &output0, vec3(255.f, 0.f, 0.f));
        }
        else
        {
            draw_point(curr_character_bone_anim_positions(Bone_RightUpLeg), camera, &output0, vec3(0.f, 255.f, 0.f));
            draw_point(curr_character_bone_anim_positions(Bone_RightFoot), camera, &output0, vec3(0.f, 255.f, 0.f));
            draw_point(curr_character_bone_anim_positions(Bone_RightToe), camera, &output0, vec3(0.f, 255.f, 0.f));
            // draw_point(contact_point, camera, &output0, vec3(0.f, 0.f, 255.f));
        }
        draw_point(curr_character_bone_anim_positions(Bone_LeftShoulder), camera, &output0, vec3(0.f, 0.f, 255.f));
        draw_point(curr_character_bone_anim_positions(Bone_LeftArm), camera, &output0, vec3(0.f, 0.f, 255.f));
        draw_point(curr_character_bone_anim_positions(Bone_LeftForeArm), camera, &output0, vec3(0.f, 0.f, 255.f));
        draw_point(curr_character_bone_anim_positions(Bone_LeftHand), camera, &output0, vec3(0.f, 0.f, 255.f));

        cv::imshow("MOOLIB", fbo_to_img(&output0));

        // std::cout << "Frame: " << t
        //           << "  vel: " << vel.x << ", " << vel.y << ", " << vel.z 
        //           << "  avel: " << avel.x << ", " << avel.y << ", " << avel.z 
        //           << "  ltoe: " << ltoe.x << ", " << ltoe.y << ", " << ltoe.z 
        //           << "  rtoe: " << rtoe.x << ", " << rtoe.y << ", " << rtoe.z 
        //           << "  contact: " << contact
        //           << std::endl;

        features(t, 0) = vel;
        features(t, 1) = avel;
        features(t, 2) = ltoe;
        features(t, 3) = rtoe;

        std::cout << " Time: " << t << " Frame: " << frame_num 
        //           << " dir: " << camera.direction.x << ", " << camera.direction.y << ", " << camera.direction.z
        //           << "  vel: " << vel.x << ", " << vel.y << ", " << vel.z 
        //           << "  avel: " << avel.x << ", " << avel.y << ", " << avel.z 
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

    // for(int t = 0; t < motion.nframes(); t++)
    // {
    //     outfile << features(t, 0).x << " " << features(t, 0).y << " " << features(t, 0).z << " " 
    //             << features(t, 1).x << " " << features(t, 1).y << " " << features(t, 1).z << " " 
    //             << features(t, 2).x << " " << features(t, 2).y << " " << features(t, 2).z << " "
    //             << features(t, 3).x << " " << features(t, 3).y << " " << features(t, 3).z << "\n";
    // }

    // outfile.close();

    // std::ifstream infile;
    // infile.open("../resources/output.txt");

    // for(int i = 0; i < motion.nframes(); i++)
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
    // db.velocities = velocities;
    // db.features = features;

    // std::cout << frames.size << ", " << velocities.size << ", " << features.rows << ", " << features.cols << std::endl;

    // // Database_sort(db, features);

    // std::cout << db.frames.size << ", " << db.velocities.size << ", " << db.features.rows << ", " << db.features.cols << std::endl;

    // Database_save(db, "../resources/database.bin");

    return 0;
}
