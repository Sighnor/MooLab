#include "motion.hpp"
#include "character.hpp"
#include "controller.hpp"
#include "database.hpp"
#include "motion.hpp"
#include "nnet.hpp"
#include "renderer.hpp"
#include "robot.hpp"

int main(int argc, char** argv)
{
    //相机
    Camera camera(1, 0.1f, 0.1f, 0.05f, 50.f, vec3(0.f, 1.2f, - 3.f), mat3());
    camera.set_view(vec3(0.f, 0.f, 1.f));
    camera.update_orientation(vec3(0.f, 1.f, 0.f));
    //控制器
    Controller controller;
    bind_controller(controller, &camera.position, &camera.direction);
    //角色，保存原始所有信息
    Character character;
    character_load(character, "../resources/character.bin");
    //机器人
    Robot robot;
    robot_load(robot, "../resources/robot.bin");
    //对角色数据的拷贝
    Mesh mesh1 = make_character_rest_mesh(character);
    Mesh mesh2 = make_robot_rest_mesh(robot);
    //材质
    Material material0(PBR);
    material0.roughness = 0.1f;
    material0.BRDFLut = img_to_tex("../resources/BRDFLut.png");
    material0.EavgLut = img_to_tex("../resources/EavgLut.png");
    Material material1(PBR);
    material1.roughness = 0.8f;
    material1.BRDFLut = img_to_tex("../resources/BRDFLut.png");
    material1.EavgLut = img_to_tex("../resources/EavgLut.png");
    //model存储mesh的拷贝，但由于mesh本身存储的是地址，实际上仍为mesh地址上的数据
    Model model1 = mesh_material_to_model(mesh1, material0);
    Model model2 = mesh_material_to_model(mesh2, material1);
    //数据库
    Database database;
    Database_load(database, "../resources/my_database.bin");
    //动作库
    BVH_Motion motion, motion1, motion2, motion3, motion4, motion5, motion6;
    Motion_load(motion, "../resources/long_motion.bin");
    motion1 = motion_sub_sequence(motion, 0, 50);
    motion2 = motion_sub_sequence(motion, 300, 350);
    motion3 = translation_and_rotation(motion1, 0, vec3(0.f, 0.5f, 0.f), vec3(0.f, 0.f, -1.f));
    motion4 = translation_and_rotation(motion2, 0, vec3(0.f, 0.5f, 0.f), vec3(-1.f, 0.f, 0.f));
    motion5 = motion_concatenate(motion3, motion4);
    //点光源
    Point_Light light(vec3(0.f, 1.f, 1.5f), vec3(0.f, 0.f, -1.f), vec3(20.f));
    Model model0 = get_light_model(light);
    //更新信息
    updated_paramters par;
    par.light_dir = &light.light_direction;
    par.light_pos = &light.light_position;
    par.light_radiance = &light.light_radiance;
    par.shadow_map = light.shadow_map;
    //渲染器
    Renderer renderer;

    renderer.models.push_back(&model0);
    renderer.models.push_back(&model1);
    renderer.models.push_back(&model2);
    renderer.point_lights.push_back(&light);

    FBO output0(540, 540);
    FBO output1(540, 540);

    //神经网络
    nnet compressor, decompressor, rot_quat, quat_rot;
    nnet_evaluation compressor_evaluation, decompressor_evaluation, rot_quat_evaluation, quat_rot_evaluation;
    nnet_load(compressor, "../resources/compressor.bin");
    nnet_load(decompressor, "../resources/decompressor.bin");
    nnet_load(rot_quat, "../resources/rot_quat.bin");
    nnet_load(quat_rot, "../resources/quat_rot.bin");
    compressor_evaluation.resize(compressor);
    decompressor_evaluation.resize(decompressor);
    rot_quat_evaluation.resize(rot_quat);
    quat_rot_evaluation.resize(quat_rot);

    // std::cout << compressor_evaluation.layers[0].size << ", " << compressor_evaluation.layers[compressor.weights.size()].size << std::endl;
    // std::cout << decompressor_evaluation.layers[0].size << ", " << decompressor_evaluation.layers[decompressor.weights.size()].size << std::endl;
    // std::cout << compressor_evaluation.layers.size() << ", " << decompressor_evaluation.layers.size() << std::endl;

    int key = 0;
    int t = 0;

    array1d<vec3> character_bone_anim_positions(motion.nbones());
    array1d<quat> character_bone_anim_rotations(motion.nbones());

    array1d<vec3> robot_bone_anim_positions1(robot.nbones());
    array1d<quat> robot_bone_anim_rotations1(robot.nbones());
    array1d<vec3> robot_bone_anim_positions2(robot.nbones());
    array1d<quat> robot_bone_anim_rotations2(robot.nbones());

    array2d<float> rotations;
    array2d<float> positions;
    array2d<quat> quaternions;
    array1d<quat> test;

    FILE* fr = fopen("../resources/rotations.bin", "rb");
    FILE* fp = fopen("../resources/positions.bin", "rb");
    FILE* fq = fopen("../resources/quaternions.bin", "rb");
    FILE* ft = fopen("../resources/test.bin", "rb");
    assert(fr != NULL && fp != NULL && fq != NULL && ft != NULL);

    array2d_read(rotations, fr);
    array2d_read(positions, fp);
    array2d_read(quaternions, fq);
    array1d_read(test, ft);

    fclose(fr);
    fclose(fp);
    fclose(fq);
    fclose(ft);

    // array1d<float> rotations(4);
    // array1d<float> positions(8);

    // std::vector<float> vals(2500);

    while (key != 27)
    {
        if(t == motion.nframes())
        {
            t = 0;
        }

        output0.set(vec3(100.f, 100.f, 200.f), 100.f);
        output1.set(vec3(100.f, 100.f, 200.f), 100.f);

        // output0.set(vec3(0.f, 0.f, 0.f), 100.f);
        // output1.set(vec3(0.f, 0.f, 0.f), 100.f);

        // camera.position.x = sin(deg_to_rad(t));
        // camera.position.y = sin(deg_to_rad(2 * t + 45));
        // camera.position.z = 2.f + 0.5f * sin(deg_to_rad(t));

        batch_forward_kinematics(motion, t, character_bone_anim_positions, character_bone_anim_rotations);

        deform_character_anim_mesh(character, character_bone_anim_positions, character_bone_anim_rotations, mesh1);

        // light.light_position = character_bone_anim_positions(0) + vec3(0.5 * sin(deg_to_rad(t)), 1.f + 0.5 * sin(deg_to_rad(2 * t)), 1.f * sin(deg_to_rad(3 * t)));
        // light.light_direction = vec3(0.f, 1.2f, 0.f) - light.light_position;
        light.light_position = vec3(5.f, 5.f, 0.f);
        renderer.models[0]->transform = mat4(eye3(), renderer.point_lights[0]->light_position) * 
                                        mat4(0.01f);
        // renderer.models[1]->transform = mat4(2.f);

        // renderer.models[2]->transform = mat4(mat3(vec3(0.5f, 0.f, 0.f), vec3(0.f, 1.f, 0.f), vec3(0.f, 0.f, 0.5f)), vec3(0.f, 0.f, 0.f));

        // renderer.models[2]->transform = mat4(Rodrigues(2 * t, vec3(0.f, -1.f, 0.f)), vec3());
        // renderer.models[2]->transform = mat4(eye3(), vec3(0.f, 0.f, 0.f));

        // slice1d_set(
        //     positions, 
        //     std::vector<float>{
        //         0.0f * sin(deg_to_rad(1.0f * t)), 
        //         0.6f * cos(deg_to_rad(1.2f * t)), 
        //         0.0f * sin(deg_to_rad(1.7f * t)), 
        //         0.5f * cos(deg_to_rad(2.1f * t)), 
        //         0.0f * cos(deg_to_rad(3.1f * t)), 
        //         0.0f * sin(deg_to_rad(4.2f * t)), 
        //         1.2f * sin(deg_to_rad(2.3f * t)),
        //         PI * cos(deg_to_rad(0.35f * t))});

        // slice1d_set(
        //     positions(t), 
        //     std::vector<float>{
        //         0.037476, 
        //         0.156693, 
        //         0.187904, 
        //         0.132341, 
        //         0.0634457, 
        //         0.0843131, 
        //         -0.247184, 
        //         3.05972});
        
        slice1d_set(
            rotations(t), 
            std::vector<float>{
                0.6f * cos(deg_to_rad(1.2f * t)), 
                0.5f * cos(deg_to_rad(2.1f * t)), 
                1.2f * sin(deg_to_rad(2.3f * t)),
                PI * cos(deg_to_rad(0.35f * t))});

        // slice1d_set(
        //     quaternions(t), 
        //     std::vector<quat>{
        //         quat(0.995357, 0.091460, -0.000000, 0.029991), 
        //         quat(0.995357, 0.091460, -0.000000, 0.029991), 
        //         quat(0.995357, 0.091460, -0.000000, 0.029991), 
        //         quat(0.995357, 0.091460, -0.000000, 0.029991), 
        //         quat(0.995357, 0.091460, -0.000000, 0.029991), 
        //         quat(0.995357, 0.091460, -0.000000, 0.029991), 
        //         quat(0.995357, 0.091460, -0.000000, 0.029991)});

        // quaternions(t, 0) = quat(1.f, 0.f, 0.f, 0.f);
        // quaternions(t, 1) = quat(-40.f, vec3(1.f, 0.f, 0.f));
        // quaternions(t, 2) = quat(1.f, 0.f, 0.f, 0.f);
        // quaternions(t, 3) = quat(1.f, 0.f, 0.f, 0.f);
        // quaternions(t, 4) = quat(1.f, 0.f, 0.f, 0.f);
        // quaternions(t, 5) = quat(45.f, vec3(1.f, 0.f, 0.f));
        // quaternions(t, 6) = quat(1.f, 0.f, 0.f, 0.f);

        // quaternions(t, 0) = quat(rad_to_deg(PI / 6.f * cos(deg_to_rad(0.35f * t))), vec3(0.f, 1.f, 0.f)) * quat(1.f, 0.f, 0.f, 0.f);
        // quaternions(t, 1) = quat(rad_to_deg(PI / 6.f * cos(deg_to_rad(0.35f * t))), vec3(0.f, 1.f, 0.f)) * quat(rad_to_deg(0.6f * cos(deg_to_rad(1.2f * t))), vec3(1.f, 0.f, 0.f));
        // quaternions(t, 2) = quat(rad_to_deg(PI / 6.f * cos(deg_to_rad(0.35f * t))), vec3(0.f, 1.f, 0.f)) * quat(1.f, 0.f, 0.f, 0.f);
        // quaternions(t, 3) = quat(rad_to_deg(PI / 6.f * cos(deg_to_rad(0.35f * t))), vec3(0.f, 1.f, 0.f)) * quat(rad_to_deg(0.5f * cos(deg_to_rad(2.1f * t))), vec3(1.f, 0.f, 0.f));
        // quaternions(t, 4) = quat(rad_to_deg(PI / 6.f * cos(deg_to_rad(0.35f * t))), vec3(0.f, 1.f, 0.f)) * quat(1.f, 0.f, 0.f, 0.f);
        // quaternions(t, 5) = quat(rad_to_deg(PI / 6.f * cos(deg_to_rad(0.35f * t))), vec3(0.f, 1.f, 0.f)) * quat(1.f, 0.f, 0.f, 0.f);
        // quaternions(t, 6) = quat(rad_to_deg(PI / 6.f * cos(deg_to_rad(0.35f * t))), vec3(0.f, 1.f, 0.f)) * quat(rad_to_deg(1.2f * sin(deg_to_rad(2.3f * t))), vec3(1.f, 0.f, 0.f));

        // slice1d_set(rotations, 
        //         std::vector<float>{
        //             cos(deg_to_rad(2.2 * t)), 
        //             circulate_float(0.1 * t, -180.f, 180.f) / 180.f,
        //             0.f});

        // slice1d_set(
        //     positions, 
        //     std::vector<float>{-0.228821, 0.322403, 0.168577, 0.09206, 0.1832901, 0.1411806, -0.109454, 2.99235});

        // nnet_evaluate(
        //     rotations(t),
        //     positions(t),
        //     compressor_evaluation,
        //     compressor);

        // nnet_evaluate(
        //     positions(t),
        //     rotations(t),
        //     decompressor_evaluation,
        //     decompressor);

        // nnet_evaluate(
        //     quaternions(t),
        //     rotations(t),
        //     quat_rot_evaluation,
        //     quat_rot);

        // nnet_evaluate(
        //     rotations(t),
        //     quaternions(t),
        //     rot_quat_evaluation,
        //     rot_quat);

        // rotations(t, 2) = rotations(t, 2) * 0.5f;

        robot_evaluate(
            rotations(t),
            positions(t));

        robot_forward_kinematics_positions(robot, robot_bone_anim_positions1, robot_bone_anim_rotations1, positions(t));
        // robot_forward_kinematics_quaternions(robot, robot_bone_anim_positions1, robot_bone_anim_rotations1, quaternions(t));
        deform_robot_anim_mesh(robot, robot_bone_anim_positions1, robot_bone_anim_rotations1, mesh2);

        // camera.set_pos(vec3(2.f, 3.f, 0.f));
        // camera.set_pos(vec3(4.f, 2.f, -0.5f));
        camera.set_pos(vec3(8.f, 2.f, 0.f));
        camera.set_view(vec3(-1.f, 0.f, 0.f));
        camera.update_orientation(vec3(0.f, 1.f, 0.f));

        // camera.set_pos(vec3(2.5f, 2.f, 0.8f));
        // camera.set_view(vec3(-1.f, 0.f, 0.f));
        // camera.update_orientation(vec3(0.f, 1.f, 0.f));

        // controller.pos_pid_control(1.f, 0.01f, 0.1f, 200.f, 1.f / 60.f, character_bone_anim_positions(0) + character_bone_anim_rotations(0) * vec3(0.f, 1.2f, -2.f));
        // controller.dir_pid_control(0.7, 0.01f, 0.5f, 200.f, 1.f / 60.f, character_bone_anim_rotations(0) * vec3(0.f, 0.f, 1.f));
        // camera.update_orientation(vec3(0.f, 1.f, 0.f));

        // vec3 yxz = quat_to_euler_YXZ(character_bone_anim_rotations(0));
        // mat3 R = Rodrigues(yxz.x, vec3(0.f, 1.f, 0.f));
        // camera.set_pos(character_bone_anim_positions(0) + character_bone_anim_rotations(0) * vec3(0.f, 1.2f, -2.f));
        // camera.set_view((character_bone_anim_rotations(0) * vec3(0.f, 0.f, 1.f)));
        // camera.update_orientation(vec3(0.f, 1.f, 0.f));
        // camera.set_view(R * vec3(0.f, 0.f, -1.f), vec3(0.f, 1.f, 0.f));

        renderer.models[1]->transform = mat4(eye3(), vec3(0.f, 0.f, -2.f)) * mat4(1.2f);
        draw(renderer.models[1], camera, &output0, &par, true);

        renderer.models[1]->transform = mat4(eye3(), vec3(0.f, 0.f, -4.f)) * mat4(1.1f);
        draw(renderer.models[1], camera, &output0, &par, true);

        renderer.models[1]->transform = mat4(eye3(), vec3(0.f, 0.f, 2.f)) * mat4(1.0f);
        draw(renderer.models[1], camera, &output0, &par, true);

        renderer.models[1]->transform = mat4(eye3(), vec3(0.f, 0.f, 4.f)) * mat4(0.9f);
        draw(renderer.models[1], camera, &output0, &par, true);

        renderer.models[1]->transform = mat4(eye3(), vec3(0.f, 0.f, 0.f)) * mat4(0.8);
        draw(renderer.models[0], camera, &output0, &par, true);
        draw(renderer.models[1], camera, &output0, &par, true);
        // draw(renderer.models[2], camera, &output0, &par, true);

        renderer.models[1]->transform = mat4(1.f);

        slice1d_set(
            rotations(t), 
            std::vector<float>{
                0.6f * cos(deg_to_rad(1.2f * t)), 
                0.5f * cos(deg_to_rad(2.1f * t)), 
                1.2f * sin(deg_to_rad(2.3f * t)),
                PI * cos(deg_to_rad(0.35f * t))});


        // slice1d_set(rotations, 
        //         std::vector<float>{
        //             cos(deg_to_rad(2.2 * t)), 
        //             circulate_float(0.1 * t, -180.f, 180.f) / 180.f,
        //             0.f});

        // slice1d_set_positions(positions, std::vector<float>(7, deg_to_rad(12.f * cos(deg_to_rad(2.2 * t)))), 2.99235);

        // nnet_evaluate(
        //     positions,
        //     rotations,
        //     decompressor_evaluation,
        //     decompressor);

        // nnet_evaluate(
        //     rotations(t),
        //     positions(t),
        //     compressor_evaluation,
        //     compressor);

        // nnet_evaluate(
        //     rotations(t),
        //     quaternions(t),
        //     rot_quat_evaluation,
        //     rot_quat);

        // nnet_evaluate(
        //     quaternions(t),
        //     rotations(t),
        //     quat_rot_evaluation,
        //     quat_rot);

        robot_evaluate(
            rotations(t),
            positions(t));

        // quaternions(t, 6) = quaternions(t, 6) * quat(rad_to_deg(-0.6f * sin(deg_to_rad(2.3f * t))), vec3(1.f, 0.f, 0.f));

        robot_forward_kinematics_positions(robot, robot_bone_anim_positions2, robot_bone_anim_rotations2, positions(t));
        // robot_forward_kinematics_quaternions(robot, robot_bone_anim_positions2, robot_bone_anim_rotations2, quaternions(t));
        deform_robot_anim_mesh(robot, robot_bone_anim_positions2, robot_bone_anim_rotations2, mesh2);

        // camera.set_pos(vec3(2.f, 3.f, -1.8f));
        // camera.set_pos(vec3(4.f, 2.f, -2.5f));
        camera.set_pos(vec3(8.f, 2.f, 0.f));
        camera.set_view(vec3(-1.f, 0.f, 0.f));
        camera.update_orientation(vec3(0.f, 1.f, 0.f));

        // camera.set_pos(vec3(2.5f, 2.f, -0.2f));
        // camera.set_view(vec3(-1.f, 0.f, 0.f));
        // camera.update_orientation(vec3(0.f, 1.f, 0.f));

        draw(renderer.models[0], camera, &output1, &par, true);
        draw(renderer.models[1], camera, &output1, &par, true);
        // draw(renderer.models[2], camera, &output1, &par, true);

        cv::imshow("left", fbo_to_img(&output0));
        cv::imshow("right", fbo_to_img(&output1));

        // char left_name[50];
        // char right_name[50];
        // sprintf(left_name, "../../Shape_Feature_Extraction/resources/img/zleft%d.png", t);
        // sprintf(right_name, "../../Shape_Feature_Extraction/resources/img/zright%d.png", t);
        // cv::imwrite(left_name, fbo_to_img(&output0));
        // cv::imwrite(right_name, fbo_to_img(&output1));

        // vals[t] = compute_val(robot_bone_anim_positions1, robot_bone_anim_positions2);

        // std::cout << vals[t] << std::endl;

        std::cout << t << std::endl;
        t++;
        // print(character_bone_anim_positions(0));
        key = cv::waitKey(1);
    }

    // std::ofstream outfile;
    // outfile.open("../resources/val_control.txt");

    // for(int i = 0; i < 2500; i++)
    // {
    //     outfile << vals[i] << "\n";
    // }

    // outfile.close();

    // for(int t = 0; t < 90; t++)
    // {
    //     // quat q = quat(t, vec3(sin(t), cos(t / 2), tan(t / 3)));

    //     // print(inv_quat(q) * q);

    //     // vec3 v1 = normalize(vec3(t / 2.f, sin(deg_to_rad(t)) / 3.f, t / 4.f));
    //     // vec3 v2 = normalize(vec3(t / 2.f, cos(deg_to_rad(t)) / 3.f, t / 4.f));
    //     // print(quat_to_euler_YXZ(euler_YXZ_to_quat(v.x, v.y, v.z)) - v);

    //     // quat q1 = quat(t, vec3(sin(deg_to_rad(t)), cos(deg_to_rad(t)), 1.f));
    //     // quat q2 = quat(2 * t, vec3(1.f , sin(deg_to_rad(t)), 1.f));

    //     // mat3 R1 = Rodrigues(t, vec3(sin(deg_to_rad(t)), cos(deg_to_rad(t)), 1.f));
    //     // mat3 R2 = Rodrigues(2 * t, vec3(1.f , sin(deg_to_rad(t)), 1.f));

    //     // std::cout << dot(q1 * v1, q2 * v1) - dot(q1, q2) << std::endl; 

    //     // print(angle_to_dir(dir_to_angle(v1)) - v1);

    //     std::cout << atan2(sin(deg_to_rad(4 * t)), cos(deg_to_rad(4 * t))) << std::endl;

    //     cv::waitKey(100);
    // }
    
    return 0;
}
