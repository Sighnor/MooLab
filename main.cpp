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
    Camera camera(1, 0.1f, 0.1f, 0.05f, 50.f, vec3(0.f, 1.3f, 5.f), mat3());
    camera.set_view(vec3(0.f, 0.f, -1.f));
    camera.update_orientation(vec3(0.f, 1.f, 0.f));
    //控制器
    Controller controller;
    bind_controller(controller, &camera.position, &camera.direction);
    //角色，保存原始所有信息
    Character character;
    character_load(character, "../resources/character.bin");
    //对角色数据的拷贝
    Mesh mesh1 = make_character_rest_mesh(character);
    //材质
    Material material0(PBR);
    material0.roughness = 0.1f;
    material0.BRDFLut = img_to_tex("../resources/BRDFLut.png");
    material0.EavgLut = img_to_tex("../resources/EavgLut.png");
    //model存储mesh的拷贝，但由于mesh本身存储的是地址，实际上仍为mesh地址上的数据
    Model model1 = mesh_material_to_model(mesh1, material0);
    //动作库
    BVH_Motion motion, motion1, motion2, motion3, motion4, motion5, motion6;
    Motion_load(motion, "../resources/long_motion.bin");
    motion1 = motion_sub_sequence(motion, 0, 250);
    motion2 = motion_sub_sequence(motion, 700, 1000);
    motion3 = translation_and_rotation(motion1, 0, vec3(0.f, 0.f, 0.f), vec3(0.f, 0.f, -1.f));
    motion4 = translation_and_rotation(motion2, 0, vec3(0.f, 0.f, 0.f), vec3(-1.f, 0.f, 0.f));
    motion5 = concatenate_two_motions(motion3, motion4, 150, 100);
    motion6 = motion;
    // motion6 = motion_sub_sequence(motion, 0, 281);
    // motion6 = motion_sub_sequence(motion, 281, 562);
    // motion6 = motion_sub_sequence(motion, 562, 13153);
    // motion6 = motion_sub_sequence(motion, 13153, 25744);
    // motion6 = motion_sub_sequence(motion, 25744, 39622);
    // motion6 = motion_sub_sequence(motion, 39622, 53500);
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
    renderer.point_lights.push_back(&light);

    FBO output0(540, 540);

    int key = 0;
    int t = 16000;

    array1d<vec3> curr_character_bone_anim_positions(motion.nbones());
    array1d<quat> curr_character_bone_anim_rotations(motion.nbones());

    array1d<vec3> last_character_bone_anim_positions(motion.nbones());
    array1d<quat> last_character_bone_anim_rotations(motion.nbones());
    // vel avel ltoe rtoe
    array2d<vec3> features(motion.nframes(), 4);
    float contact = 0;

    while (key != 27)
    {
        if(t == motion.nframes())
        {
            t = 0;
            // BVH_Motion motion = translation_and_rotation(motion, 0, curr_character_bone_anim_positions(Bone_Entity), curr_character_bone_anim_rotations(Bone_Entity) * vec3(0.f, 0.f, 1.f));
        }

        // motion5 = translation_and_rotation(motion5, 0, curr_character_bone_anim_positions(Bone_Entity), curr_character_bone_anim_rotations(Bone_Entity) * vec3(0.f, 0.f, 1.f));

        output0.set(vec3(100.f, 100.f, 200.f), 100.f);

        last_character_bone_anim_positions = curr_character_bone_anim_positions;
        last_character_bone_anim_rotations = curr_character_bone_anim_rotations;

        batch_forward_kinematics_full(
            motion, 
            t, 
            curr_character_bone_anim_positions, 
            curr_character_bone_anim_rotations);

        // batch_forward_kinematics_part(
        //     motion, 
        //     t, 
        //     curr_character_bone_anim_positions, 
        //     curr_character_bone_anim_rotations, 
        //     vec3(0.f), 
        //     quat(0.f, vec3(0.f, 0.f, 1.f)));

        deform_character_anim_mesh(
            character, 
            curr_character_bone_anim_positions, 
            curr_character_bone_anim_rotations, 
            mesh1);

        features(t, 0) = inv_quat(curr_character_bone_anim_rotations(Bone_Entity)) 
                   * vec_to_vel(last_character_bone_anim_positions(Bone_Entity), curr_character_bone_anim_positions(Bone_Entity));

        features(t, 1) = inv_quat(curr_character_bone_anim_rotations(Bone_Entity)) 
                    * quat_to_avel(last_character_bone_anim_rotations(Bone_Entity), curr_character_bone_anim_rotations(Bone_Entity));

        features(t, 2) = inv_quat(curr_character_bone_anim_rotations(Bone_Entity)) 
                    * (curr_character_bone_anim_positions(Bone_LeftToe) - curr_character_bone_anim_positions(Bone_Entity));

        features(t, 3) = inv_quat(curr_character_bone_anim_rotations(Bone_Entity)) 
                    * (curr_character_bone_anim_positions(Bone_RightToe) - curr_character_bone_anim_positions(Bone_Entity));

        if(t > 0)
        {
            vec3 lvel = vec_to_vel(features(t - 1, 2), features(t, 2));
            vec3 rvel = vec_to_vel(features(t - 1, 3), features(t, 3));
            if(dot(features(t, 0), lvel) < 0.f)
            {
                contact = 0;
            }
            else
            {
                contact = 1;
            }
        }

        light.light_position = curr_character_bone_anim_positions(Bone_Entity) + 
                               curr_character_bone_anim_rotations(Bone_Entity) * vec3(0.f, 1.5f, -1.5f);
        light.light_radiance = 20.f + 10.f * sin(deg_to_rad(0.5f * t));

        controller.pos_pid_control(
            1.f, 
            0.01f, 
            0.1f, 
            500.f, 
            1.f / 60.f, 
            curr_character_bone_anim_positions(Bone_Entity) + curr_character_bone_anim_rotations(Bone_Entity) * vec3(0.f, 1.2f, -2.f));
        controller.dir_pid_control(
            0.7, 
            0.01f, 
            0.1f, 
            500.f, 
            1.f / 60.f, 
            curr_character_bone_anim_rotations(Bone_Entity) * vec3(0.f, 0.f, 1.f));
        camera.update_orientation(vec3(0.f, 1.f, 0.f));

        renderer.models[0]->transform = mat4(eye3(), renderer.point_lights[0]->light_position) * mat4(0.01f);
        // renderer.models[0]->transform = mat4(eye3(), curr_character_bone_anim_positions(Bone_Entity)) * mat4(0.1f);

        draw(renderer.models[0], camera, &output0, &par, true);
        draw(renderer.models[1], camera, &output0, &par, true);
        if(contact == 0)
        {
            draw_point(curr_character_bone_anim_positions(Bone_LeftToe), camera, &output0, vec3(255.f, 0.f, 0.f));
        }
        else
        {
            draw_point(curr_character_bone_anim_positions(Bone_RightToe), camera, &output0, vec3(255.f, 0.f, 0.f));
        }

        cv::imshow("MOOLIB", fbo_to_img(&output0));

        std::cout << "Frame: " << t
                  << "    vel: " << features(t, 0).x << ", " << features(t, 0).y << ", " << features(t, 0).z 
                  << "    avel: " << features(t, 1).x << ", " << features(t, 1).y << ", " << features(t, 1).z 
                  << "    ltoe: " << features(t, 2).x << ", " << features(t, 2).y << ", " << features(t, 2).z 
                  << "    rtoe: " << features(t, 3).x << ", " << features(t, 3).y << ", " << features(t, 3).z 
                  << "    contact: " << contact
                  << std::endl;

        // std::cout << "Frame: " << t << std::endl;
        t++;
        key = cv::waitKey(1);
    }

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

    // FILE* fft = fopen("../resources/features.bin", "wb");
    // assert(fft != NULL);
    // array2d_write(features, fft);
    // fclose(fft);

    return 0;
}
