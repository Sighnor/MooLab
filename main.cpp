#include "character.hpp"
#include "controller.hpp"
#include "database.hpp"
#include "motion.hpp"
#include "renderer.hpp"

int main(int argc, char** argv)
{
	// mat3 m(vec3(2.f, 0.f, 0.f), vec3(0.f, 2.f, 0.f), vec3(0.f, 0.f, 5.f));
	// print(m);
	// print(inv_mat(m));

	// mat3 R = Rodrigues(vec3(0, 0, 1), 60);

	// Camera camera(1, 1920, 1080, 0.1f, 20.f, vec3(1.f, 2.f, 3.f), R);

	// print(camera.position);
	// print(camera.orientation);
	// print(mat4(camera.orientation, camera.position));
	// print(camera.get_inverse_view_matrix() * camera.get_view_matrix());
	// print(camera.get_inverse_projection_matrix() * camera.get_projection_matrix());

	// float alpha = 0.f;
	// float beta = 0.f;
	// float gamma = 0.f;
	// bool flag;

	// vec3 v = vec3(0.2f, 2.f, 0.f);
	// vec3 x = vec3(-1.f, 0.f, 0.f);
	// vec3 y = vec3(1.f, 0.f, 0.f);
	// vec3 z = vec3(0.f, 3.f, 0.f);

	// flag = inside_3dtriangle(x, y, z, v, alpha, beta, gamma);

	// printf("%d, %f, %f, %f\n", flag, alpha, beta, gamma);

	// print(normalize(cross(y - x, z - x)));

	// mat3 m = get_coordinate_matrix(vec3(2.f, 1.f, 1.f), vec3(0.f, 1.f, 0.f));

	// print(m);

	// flag = inside_2dtriangle(vec3_to_vec2(m * x), vec3_to_vec2(m * y), vec3_to_vec2(m * z), vec3_to_vec2(m * v), alpha, beta, gamma);

	// printf("%d, %f, %f, %f\n", flag, alpha, beta, gamma);

	// flag = inside_3dtriangle(m * x, m * y, m * z, m * v, alpha, beta, gamma);

	// printf("%d, %f, %f, %f\n", flag, alpha, beta, gamma);

	// Renderer renderer;

	// updated_paramters par;

	// draw(renderer.models[0], camera, renderer.point_lights[0].shadow_map, NULL);

	// draw(renderer.models[1], camera, camera.fbo, NULL);
 
	// draw(renderer.models[2], camera, NULL, &par);

	//相机
	Camera camera(1, 0.1f, 0.1f, 0.05f, 50.f, vec3(0.f, 1.2f, 1.f), mat3());
	camera.set_view(vec3(0.f, 0.f, -1.f), vec3(0.f, 1.f, 0.f));
	// print(camera.orientation);
	//角色，保存原始所有信息
	Character character;
	character_load(character, "../resources/character.bin");
	//对角色数据的拷贝
	Mesh mesh = make_character_rest_mesh(character);
	// Mesh mesh = triangle(eye3(), vec3());
	//材质
	Material material;
	material.roughness = 1.f;
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

	FBO output(540, 540, vec3(0.f, 0.f, 0.f));

	int key = 0;
	int t = 0;

	while (key != 27)
	{
		// camera.position.x = sin(deg_rad(t));
		// camera.position.y = sin(deg_rad(2 * t + 45));
		renderer.models[0]->transform = mat4(Rodrigues(vec3(0.f, 1.f, 0.f), t), vec3());
		draw(renderer.models[0], camera, &output, &par);
		cv::imshow("img", FBO_to_IMG(&output));
		t++;
		std::cout << t << std::endl;
		key = cv::waitKey(1);
		output.colors.set(vec3(0.f));
	}
	
	// draw(renderer.models[0], camera, &output, &par);
	// cv::imshow("img", FBO_to_IMG(&output));
	// t++;
	// std::cout << t << std::endl;
	// key = cv::waitKey(100000);

	
	return 0;
}
