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
	//角色，保存原始所有信息
	Character character;
	character_load(character, "../resources/character.bin");
	//对角色数据的拷贝
	Mesh mesh = make_character_rest_mesh(character);
	//材质
	Material material;
	material.roughness = 0.5f;
	//model存储mesh的拷贝，但由于mesh本身存储的是地址，实际上仍为mesh地址上的数据
	Model model = mesh_material_to_model(mesh, material);
	//点光源
	Point_Light light(vec3(0.f, 5.f, 5.f), vec3(0.f, -1.f, -1.f), vec3(20.f));
	//更新信息
	updated_paramters par;
	par.light_dir = &light.light_direction;
	par.light_pos = &light.light_position;
	par.light_radiance = &light.light_radiance;
	par.shadow_map = light.shadow_map;
	//渲染器
	Renderer renderer;

	renderer.models.push_back(&model);
	renderer.point_lights.push_back(&light);

	FBO output(540, 540, vec3(100.f, 100.f, 200.f));

	int key = 0;
	int t = 0;

	while (key != 27)
	{
		// camera.position.x = sin(deg_rad(t));
		// camera.position.y = sin(deg_rad(2 * t + 45));
		camera.position.z = 2.f + 0.5f * sin(deg_rad(t));
		renderer.models[0]->transform = mat4(Rodrigues(vec3(0.f, 1.f, 0.f), t), vec3());

		draw(renderer.models[0], camera, &output, &par);
		cv::imshow("img", FBO_to_IMG(&output));
		t++;
		std::cout << t << std::endl;
		key = cv::waitKey(1);
		output.colors.set(vec3(100.f, 100.f, 200.f));
	}
	
	return 0;
}
