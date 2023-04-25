#include "camera.hpp"
#include "controller.hpp"
#include "database.hpp"
#include "light.hpp"
#include "mesh.hpp"
#include "model.hpp"
#include "renderer.hpp"
#include "shader.hpp"
#include "matrix.hpp"

int main(int argc, char** argv)
{
	mat3 m(vec3(2.f, 0.f, 0.f), vec3(0.f, 2.f, 0.f), vec3(0.f, 0.f, 5.f));
	print(m);
	print(inv_mat(m));

	mat3 R = Rodrigues(vec3(0, 0, 1), 45);
	print(R);

	Camera camera(1, 1, 1, 0.1f, 10.f, vec3(), mat3());

	Renderer renderer;

	updated_paramters par;

	draw(renderer.models[0], camera, renderer.point_lights[0].shadow_map, NULL);

	draw(renderer.models[1], camera, camera.fbo, NULL);
 
	draw(renderer.models[2], camera, NULL, &par);
	
	return 0;
}
