#include "character.hpp"
#include "controller.hpp"
#include "database.hpp"
#include "motion.hpp"
#include "renderer.hpp"

int main(int argc, char** argv)
{
	//相机
	Camera camera(1, 0.1f, 0.1f, 0.05f, 50.f, vec3(0.f, 1.2f, 2.f), mat3());
	camera.set_view(vec3(0.f, 0.f, -1.f), vec3(0.f, 1.f, 0.f));
	//控制器
	Controller controller;
	bind_controller(controller, &camera.position, &camera.direction);
	//角色，保存原始所有信息
	Character character;
	character_load(character, "../resources/character.bin");
	//对角色数据的拷贝
	Mesh mesh = make_character_rest_mesh(character);
	//材质
	Material material(PBR);
	material.roughness = 0.35f;
	material.BRDFLut = img_to_tex("../resources/BRDFLut.png");
	material.EavgLut = img_to_tex("../resources/EavgLut.png");
	//model存储mesh的拷贝，但由于mesh本身存储的是地址，实际上仍为mesh地址上的数据
	Model model1 = mesh_material_to_model(mesh, material);
	Model model2 = mesh_material_to_model(mesh, material);
	//数据库
	Database database;
	Database_load(database, "../resources/my_database.bin");
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

	FBO output(540, 540);

	int key = 0;
	int t = 0;

	array1d<vec3> bone_anim_positions(database.nbones());
    array1d<quat> bone_anim_rotations(database.nbones());

	while (key != 27)
	{
		if(t == database.nframes())
		{
			t = 0;
		}

		output.set(vec3(100.f, 100.f, 200.f), 100.f);

		// camera.position.x = sin(deg_to_rad(t));
		// camera.position.y = sin(deg_to_rad(2 * t + 45));
		// camera.position.z = 2.f + 0.5f * sin(deg_to_rad(t));

		batch_forward_kinematics(database, t, bone_anim_positions, bone_anim_rotations);

		deform_character_anim_mesh(character, bone_anim_positions, bone_anim_rotations, mesh);

		light.light_position = bone_anim_positions(0) + vec3(0.5 * sin(deg_to_rad(t)), 1.f + 0.5 * sin(deg_to_rad(2 * t)), 1.f * sin(deg_to_rad(3 * t)));
		light.light_direction = vec3(0.f, 1.2f, 0.f) - light.light_position;
		renderer.models[0]->transform = mat4(eye3(), renderer.point_lights[0]->light_position) * 
										mat4(0.01f);
		// renderer.models[1]->transform = mat4(Rodrigues(2 * t, vec3(0.f, -1.f, 0.f)), vec3());

		// camera.set_pos(bone_anim_positions(0) + vec3(0.f, 1.2f, 2.f));

		controller.pos_pid_control(0.5f, 0.01f, 0.1f, 1000.f, 1.f / 60.f, bone_anim_positions(0) + vec3(0.f, 1.2f, 2.f));

		// vec3 yxz = quat_to_euler_YXZ(bone_anim_rotations(0));
		// mat3 R = Rodrigues(yxz.x, vec3(0.f, 1.f, 0.f));
		// camera.set_pos(bone_anim_positions(0) - bone_anim_rotations(0) * vec3(0.f, -1.2f, 2.f));
		// camera.set_view(-(bone_anim_rotations(0) * vec3(0.f, 0.f, -1.f)), vec3(0.f, 1.f, 0.f));
		// camera.set_view(R * vec3(0.f, 0.f, -1.f), vec3(0.f, 1.f, 0.f));

		draw(renderer.models[0], camera, &output, &par, true);
		draw(renderer.models[1], camera, &output, &par, true);
		// draw(renderer.models[1], camera, &output, &par, true);
		cv::imshow("img", fbo_to_img(&output));
		t++;
		std::cout << t << std::endl;
		key = cv::waitKey(1);
	}

	// for(int t = 0; t < 360; t++)
	// {
	// 	// quat q = quat(t, vec3(sin(t), cos(t / 2), tan(t / 3)));

	// 	// print(inv_quat(q) * q);

	// 	vec3 v = vec3(t / 2.f, t / 3.f, t / 4.f);
	// 	print(quat_to_euler_YXZ(euler_YXZ_to_quat(v.x, v.y, v.z)) - v);

	// 	cv::waitKey(100);
	// }
	
	return 0;
}
