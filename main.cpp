#include "motion.hpp"
#include "character.hpp"
#include "controller.hpp"
#include "database.hpp"
#include "motion.hpp"
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
    material1.roughness = 0.5f;
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

    int key = 0;
    int t = 0;

<<<<<<< HEAD
	array1d<vec3> bone_anim_positions(motion.nbones());
        array1d<quat> bone_anim_rotations(motion.nbones());
=======
    array1d<vec3> character_bone_anim_positions(motion.nbones());
    array1d<quat> character_bone_anim_rotations(motion.nbones());
>>>>>>> 030d30a (v1.4)

    array1d<vec3> robot_bone_anim_positions(robot.nbones());
    array1d<quat> robot_bone_anim_rotations(robot.nbones());

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

        // batch_forward_kinematics(motion, t, character_bone_anim_positions, character_bone_anim_rotations);

        // deform_character_anim_mesh(character, character_bone_anim_positions, character_bone_anim_rotations, mesh1);

        robot_forward_kinematics(robot, robot_bone_anim_positions, robot_bone_anim_rotations, 12.f * cos(deg_to_rad(4.4 * t)), 12.f * cos(deg_to_rad(1.f * t)), 0.f);

        deform_robot_anim_mesh(robot, robot_bone_anim_positions, robot_bone_anim_rotations, mesh2);

        light.light_position = character_bone_anim_positions(0) + vec3(0.5 * sin(deg_to_rad(t)), 1.f + 0.5 * sin(deg_to_rad(2 * t)), 1.f * sin(deg_to_rad(3 * t)));
        light.light_direction = vec3(0.f, 1.2f, 0.f) - light.light_position;
        renderer.models[0]->transform = mat4(eye3(), renderer.point_lights[0]->light_position) * 
                                        mat4(0.01f);

        // renderer.models[2]->transform = mat4(Rodrigues(2 * t, vec3(0.f, -1.f, 0.f)), vec3());
        // renderer.models[2]->transform = mat4(eye3(), vec3(0.f, 0.f, 0.f));

        camera.set_pos(vec3(5.f, 2.f, 0.f));
        camera.set_view(vec3(-1.f, 0.f, 0.f));
        camera.update_orientation(vec3(0.f, 1.f, 0.f));

        // controller.pos_pid_control(1.f, 0.01f, 0.1f, 200.f, 1.f / 60.f, character_bone_anim_positions(0) + character_bone_anim_rotations(0) * vec3(0.f, 1.2f, -2.f));
        // controller.dir_pid_control(0.7, 0.01f, 0.5f, 200.f, 1.f / 60.f, character_bone_anim_rotations(0) * vec3(0.f, 0.f, 1.f));
        // camera.update_orientation(vec3(0.f, 1.f, 0.f));

        // vec3 yxz = quat_to_euler_YXZ(character_bone_anim_rotations(0));
        // mat3 R = Rodrigues(yxz.x, vec3(0.f, 1.f, 0.f));
        // camera.set_pos(character_bone_anim_positions(0) + character_bone_anim_rotations(0) * vec3(0.f, 1.2f, -2.f));
        // camera.set_view((character_bone_anim_rotations(0) * vec3(0.f, 0.f, 1.f)));
        // camera.update_orientation(vec3(0.f, 1.f, 0.f));
        // camera.set_view(R * vec3(0.f, 0.f, -1.f), vec3(0.f, 1.f, 0.f));

        // draw(renderer.models[0], camera, &output0, &par, true);
        // draw(renderer.models[1], camera, &output0, &par, true);

        draw(renderer.models[2], camera, &output0, &par, true);

        // camera.set_pos(vec3(5.f, 2.f, -0.5f));
        // camera.set_view(vec3(-1.f, 0.f, 0.f));
        // camera.update_orientation(vec3(0.f, 1.f, 0.f));

        robot_forward_kinematics(robot,
                                robot_bone_anim_positions,
                                robot_bone_anim_rotations,
                                12.f * cos(deg_to_rad(4.4 * t)) * (1.f + 0.15f * sin(deg_to_rad(t))),
                                12.f * cos(deg_to_rad(1.f * t)) * (1.f + 0.15f * sin(deg_to_rad(t))),
                                0.f);

        deform_robot_anim_mesh(robot, robot_bone_anim_positions, robot_bone_anim_rotations, mesh2);

        draw(renderer.models[2], camera, &output1, &par, true);

        cv::imshow("left", fbo_to_img(&output0));
        cv::imshow("right", fbo_to_img(&output1));

        // char left_name[50];
        // char right_name[50];
        // sprintf(left_name, "../resources/zleft%d.png", t);
        // sprintf(right_name, "../resources/zright%d.png", t);
        // cv::imwrite(left_name, fbo_to_img(&output0));
        // cv::imwrite(right_name, fbo_to_img(&output1));
        t++;
        std::cout << t << std::endl;
        // print(character_bone_anim_positions(0));
        key = cv::waitKey(1);
    }

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
