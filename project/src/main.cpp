/**
	Objectives:
	- Complete rigid transform interpolation in the function evaluate_local function in file skeleton.cpp
	- Complete the Linear Blend Skinning computation in the function skinning.cpp
*/

#include "vcl/vcl.hpp"
#include <iostream>


#include "forest.hpp"



using namespace vcl;

//for patronus and spell
struct particle
{
	vec3 p0;
	float t0;
	vec3 v0;
	vec3 color;
	trajectory_drawable trajectory;
};



struct gui_parameters {
	bool display_frame = true;
	
};

struct user_interaction_parameters {
	vec2 mouse_prev;
	timer_fps fps_record;
	gui_parameters gui;
	mesh_drawable global_frame;
	bool cursor_on_gui;
};
user_interaction_parameters user;

struct scene_environment
{
	camera_around_center camera;
	mat4 projection;
	vec3 light;
};
scene_environment scene;

//Forest


struct forest_structure
{
    grid_2D<vec3> position;

    grid_2D<vec3> normal;

    buffer<uint3> triangle_connectivity;
    mesh_drawable visual;
};

forest_structure forest;

//Dementor
struct dementor_structure
{
	int N_cloth;
	/*
    grid_2D<vec3> position;
    grid_2D<vec3> velocity;
    grid_2D<vec3> forces;

    grid_2D<vec3> normal;
	*/

	buffer<vec3> position;
	buffer<vec3> velocity;
	buffer<vec3> normal;
    buffer<uint3> triangle_connectivity;
    mesh_drawable visual;

    std::map<size_t,vec3> positional_constraints;

    float x0;
    float z0;
    float y0;

    float fi;
	
};

//character

mesh_drawable head;
mesh_drawable body;
mesh_drawable arm;
mesh_drawable legs;
mesh_drawable hair1;
mesh_drawable hair2;

//Patronus
mesh_drawable sphere;
mesh_drawable deer;
GLuint texture_deer;
mesh deer_mesh;

//trees
mesh_drawable tree;
mesh_drawable trunk;
mesh_drawable leaves;

vcl::buffer<vcl::vec2> positions;

//lake
mesh_drawable lake;

mesh_drawable castle;

timer_interval timer;
std::vector<particle> particle_deer_collection;
std::vector<particle> particle_spell_collection;
timer_event_periodic timer_particle(0.01f);
timer_event_periodic timer_spell(0.001f);

//dementors

std::vector<dementor_structure> dementors;

// spells
mesh_drawable sphere2;

bool levitate = false;

//interaction

void mouse_move_callback(GLFWwindow* window, double xpos, double ypos);
void window_size_callback(GLFWwindow* window, int width, int height);

void initialize_data();
void initialize_forest();
void initialize_deer();
void initialize_body();
void initialize_dementor();
void initialize_lake();

void display_scene();
void display_interface();
//void compute_deformation();
void update_new_content(mesh const& shape, GLuint texture_id);

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);


int main(int argc, char** argv)
{

	std::cout << "Run " << argv[0] << std::endl;

	int const width = 1280, height = 1024;
	GLFWwindow* window = create_window(width, height);
	window_size_callback(window, width, height);
	std::cout << opengl_info_display() << std::endl;;

	imgui_init(window);
	glfwSetCursorPosCallback(window, mouse_move_callback);
	glfwSetWindowSizeCallback(window, window_size_callback);
	glfwSetKeyCallback (window, key_callback);

	std::cout<<"Initialize data ..."<<std::endl;
	initialize_data();


	std::cout<<"Start animation loop ..."<<std::endl;
	user.fps_record.start();
	timer.start();
	timer.scale = 0.2f;
	timer_particle.start();

	scene.camera.manipulator_rotate_trackball({0,0}, {1,0});
	scene.camera.manipulator_rotate_trackball({0,0}, {1,0});
	scene.camera.manipulator_rotate_trackball({0,0}, {1,0});

	glEnable(GL_DEPTH_TEST);
	while (!glfwWindowShouldClose(window))
	{
		scene.light = scene.camera.position();
		user.fps_record.update();
		timer.update();

		glClearColor(0.5f, 0.5f, 0.5f, 0.5f);
		glClear(GL_COLOR_BUFFER_BIT);
		glClear(GL_DEPTH_BUFFER_BIT);
		imgui_create_frame();
		if(user.fps_record.event) {
			std::string const title = "VCL Display - "+str(user.fps_record.fps)+" fps";
			glfwSetWindowTitle(window, title.c_str());
		}

		ImGui::Begin("GUI",NULL,ImGuiWindowFlags_AlwaysAutoResize);
		user.cursor_on_gui = ImGui::IsAnyWindowFocused();

		if(user.gui.display_frame)
			draw(user.global_frame, scene);
		
		//glfwSetKeyCallback(window, key_callback);
		
		
		//display_interface();
		display_scene();

		ImGui::End();
		imgui_render_frame(window);
		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	imgui_cleanup();
	glfwDestroyWindow(window);
	glfwTerminate();

	return 0;
}


void initialize_data()
{
	GLuint const shader_mesh = opengl_create_shader_program(opengl_shader_preset("mesh_vertex"), opengl_shader_preset("mesh_fragment"));
	GLuint const shader_uniform_color = opengl_create_shader_program(opengl_shader_preset("single_color_vertex"), opengl_shader_preset("single_color_fragment"));
	GLuint const texture_white = opengl_texture_to_gpu(image_raw{1,1,image_color_type::rgba,{255,255,255,255}});
	mesh_drawable::default_shader = shader_mesh;
	mesh_drawable::default_texture = texture_white;
	curve_drawable::default_shader = shader_uniform_color;
	segments_drawable::default_shader = shader_uniform_color;

	user.global_frame = mesh_drawable(mesh_primitive_frame());
	user.gui.display_frame = false;
	scene.camera.distance_to_center = 5.0f;
	scene.camera.manipulator_translate_in_plane(vec2{0.0, -2.0});

	mesh new_shape;

	GLuint const texture_spell = opengl_texture_to_gpu(image_load_png("assets/sparcle.png"));
	float const L = 0.2f; //size of the quad
	sphere2 = mesh_drawable(mesh_primitive_quadrangle({-L,-L,0},{L,-L,0},{L,L,0},{-L,L,0}));
	sphere2.texture = texture_spell;

	//mountain = mesh_drawable(mesh_load_file_obj("assets/mountain.obj"));


	initialize_forest();
	initialize_body();
	initialize_deer();
	initialize_lake();


	initialize_dementor();

	generate_tree_positions(positions);
	//initialize_simulation_parameters(dementor.parameters, 1.0f, dementor.position.dimension.x);
}


void initialize_forest(){

	mesh const forest_mesh = create_terrain();
	forest.position = grid_2D<vec3>::from_buffer(forest_mesh.position, 30,30);
	forest.visual.clear();
	forest.visual = mesh_drawable(forest_mesh);
	forest.visual.texture = opengl_texture_to_gpu(image_load_png("assets/ground.png"));
	forest.visual.shading.phong = {0.3f, 0.7f, 0.05f, 32};

	forest.triangle_connectivity = forest_mesh.connectivity;
	forest.normal = grid_2D<vec3>::from_buffer(forest_mesh.normal, 30, 30);

	
	GLuint texture_trunk = opengl_texture_to_gpu(image_raw{ 1,1,image_color_type::rgba,{64,42,13,255} });
	GLuint texture_leaves = opengl_texture_to_gpu(image_load_png("assets/leaves.png"));
	tree = mesh_drawable(mesh_load_file_obj("assets/Tree.obj"));
	tree.transform.scale = 1.2;
	tree.texture = texture_leaves;

	trunk = mesh_drawable(mesh_load_file_obj("assets/trunk.obj"));
	trunk.texture = texture_trunk;



	return;
}

void initialize_lake(){

	lake = mesh_drawable(mesh_primitive_disc(20.0f, {0,-1.0f, 30.0f}, {0,1,0}, 100));
	lake.shading.color = {0.5f,0.8f,1.0f};
	lake.texture = opengl_texture_to_gpu(image_load_png("assets/lake.png"));

}

void initialize_body(){
	head = mesh_drawable(mesh_load_file_obj("assets/head.obj"));
	body = mesh_drawable(mesh_load_file_obj("assets/upper_body.obj"));
	arm = mesh_drawable(mesh_load_file_obj("assets/rightarm_w_wand.obj"));
	legs = mesh_drawable(mesh_load_file_obj("assets/legs.obj"));
	hair1 = mesh_drawable(mesh_load_file_obj("assets/hair1.obj"));
	hair2 = mesh_drawable(mesh_load_file_obj("assets/hair2.obj"));

	body.transform.translate = { 0.0f, -1.0f, 0.0f};
	body.transform.scale = 0.5;

	head.transform.translate = { 0.0f, -1.0f, 0.0f};
	head.transform.scale = 0.5;

	legs.transform.translate = { 0.0f, -1.0f, 0.0f};
	legs.transform.scale = 0.5;

	arm.transform.translate={0.2f, -0.5f, 0.0f};
	arm.transform.scale = 0.5f ;

	hair1.transform.translate = { 0.0f, -1.0f, 0.0f};
	hair1.transform.scale = 0.5;

	hair2.transform.translate = { 0.0f, -1.0f, 0.0f};
	hair2.transform.scale = 0.5;


	GLuint texture_face= opengl_texture_to_gpu(image_raw{ 1,1,image_color_type::rgba,{240, 180, 140, 255} });
	GLuint texture_body= opengl_texture_to_gpu(image_raw{ 1,1,image_color_type::rgba,{150, 69, 69, 255} });
	GLuint texture_arm= opengl_texture_to_gpu(image_raw{ 1,1,image_color_type::rgba,{150, 69, 69, 255} });
	GLuint texture_legs = opengl_texture_to_gpu(image_load_png("assets/jeans.png"));
	GLuint texture_hair = opengl_texture_to_gpu(image_raw{ 1,1,image_color_type::rgba,{0, 0, 0, 255} });

	head.texture = texture_face;
	body.texture = texture_body;
	arm.texture = texture_arm;
	legs.texture = texture_legs;
	hair1.texture = texture_hair;
	hair2.texture = texture_hair;


	

	return;
}

void initialize_deer() {

	texture_deer = opengl_texture_to_gpu(image_raw{ 1,1,image_color_type::rgba,{244,255,254,155} });

	// texture_deer = opengl_texture_to_gpu(image_load_png("assets/deer.png"));
	deer_mesh = mesh_load_file_obj("assets/deer.obj");
	int N = deer_mesh.position.size();
	rotation const r = trackball_rotation(vec2(0.0f, 0.0f), vec2(0.0f, 2.9f));
	for (int i = 0; i < N; i++) {
		deer_mesh.position[i] = 0.01 * deer_mesh.position[i];
		deer_mesh.position[i] = r * deer_mesh.position[i];
		deer_mesh.position[i] += vec3{ 0.0f, -1.0f, 5.5f };

	}

	GLuint const texture_billboard = opengl_texture_to_gpu(image_load_png("assets/sparcle.png"));
	deer = mesh_drawable(deer_mesh);
	deer.transform.scale = 1.2f;
	deer.shading.alpha = 0.5;
	// sphere = mesh_drawable(mesh_primitive_sphere(1.0f));
	deer.texture = texture_deer;

	float const L = 0.50f; //size of the quad
	sphere = mesh_drawable(mesh_primitive_quadrangle({ -L,-L,0 }, { L,-L,0 }, { L,L,0 }, { -L,L,0 }));
	sphere.texture = texture_billboard; 
	return;

}


void initialize_dementor(){
	
	for (int i=0; i<10; i++){

		float x = rand_interval(-20.0f,20.0f);
		float y = rand_interval(10.0,50.0f);
		float z = rand_interval(5.0f, 10.0f);

		float fi = rand_interval(0,3.14f);
	
		dementor_structure dementor;

		mesh const cloth_mesh = mesh_load_file_obj("assets/dementor.obj");
/*
		float const N_cloth = 30;
		mesh const cloth_mesh = mesh_primitive_grid({0,0,0},{10,0,0},{10,0,10},{0,0,10},N_cloth, N_cloth);

		

		dementor.position = grid_2D<vec3>::from_buffer(cloth_mesh.position, N_cloth, N_cloth);
		dementor.normal = grid_2D<vec3>::from_buffer(cloth_mesh.normal, N_cloth, N_cloth);

		dementor.velocity.clear();
		dementor.velocity.resize(N_cloth*N_cloth);

		dementor.forces.clear();
		dementor.forces.resize(N_cloth, N_cloth);

		dementor.positional_constraints.clear();
		dementor.positional_constraints[dementor.position.index_to_offset(0,0)] = dementor.position(0,0);
		dementor.positional_constraints[dementor.position.index_to_offset(N_cloth-1,0)] = dementor.position(int(N_cloth-1),0);
*/
		dementor.position = cloth_mesh.position;
		dementor.normal = cloth_mesh.normal;
		dementor.visual.clear();
		dementor.visual = mesh_drawable(cloth_mesh);
		dementor.visual.shading.phong = {0.3f, 0.7f, 0.05f, 32};

		dementor.triangle_connectivity = cloth_mesh.connectivity;
		
		dementor.visual.transform.scale = 0.1f;
		dementor.visual.texture = opengl_texture_to_gpu(image_load_png("assets/dementor_texture.png"));

		
		
		dementor.x0 = x;
		dementor.y0 = y;
		dementor.z0 = z;
		dementor.fi = fi;
		dementors.push_back(dementor);

	}
	

}
//spell
particle create_new_spell_particle(float t)
{
	particle particle;
	particle.t0 = t;
	particle.p0 = {-0.5f, 0.8f, 0.7f};

	float const theta = rand_interval(0,0.5*pi);
	float const fi = rand_interval(-0.2*pi,0.2*pi);
	float const v0y = rand_interval(5,10);
	particle.v0 = {v0y*std::cos(theta)*std::sin(fi),v0y*std::sin(theta), v0y*std::cos(theta)*std::cos(fi)};

	return particle;
}

//patronus
particle create_new_particle(float t)
{
	particle particle;
	particle.t0 = t;

	int rand_pos = rand() % deer_mesh.position.size();

	particle.p0 = deer_mesh.position[rand_pos];
	particle.v0 = rand_interval(0.05f, 0.3f) * deer_mesh.normal[rand_pos];
	particle.color = { 0.9f + rand_interval(0,0.05f),0.9f + rand_interval(0,0.05f), 1.0f};

	return particle;
}

vec3 compute_particle_position(particle const& particle, float t_current)
{
	float const t = t_current - particle.t0;

	return particle.p0 + particle.v0 * t;

}

vec3 compute_particle_spell_position(particle const& particle, float t_current)
{
	float const t = t_current - particle.t0;
	vec3 g_weak = {0,-5.0f,0};

	return 0.5f*g_weak*t*t + particle.v0*t + particle.p0;

}

template <typename T>
void remove_old_particles(std::vector<T>& particles, float t_current, float t_max)
{
	for (auto it = particles.begin(); it != particles.end();)
	{
		if (t_current - it->t0 > t_max)
			it = particles.erase(it);
		if (it != particles.end())
			++it;
	}
}



void draw_trees() {

	for (size_t k = 0; k< positions.size(); k++){
		float x = positions[k][0];
		float z = positions[k][1];
		tree.transform.translate = { x, -1.2f, z};
		draw(tree, scene);
		trunk.transform.translate = { x, -1.2f, z};
		draw(trunk, scene);
	}


	
}

void display_scene()
{
	//Expecto Patronum : if E is pressed


	if(levitate){
		timer_spell.update();


		//raise arm

		arm.transform.rotate = vcl::rotation({1,0,0}, -3.14f/2.5f);
		arm.transform.translate = {0.0f,0.0f,1.3f};

		//emit particles

		if (timer_spell.event)
			particle_spell_collection.push_back(create_new_spell_particle(timer_spell.t));

		timer_particle.update();
		if (timer_particle.event)
			particle_deer_collection.push_back(create_new_particle(timer_particle.t));
		
		
		deer.shading.alpha = 0.25f * (3-(timer_particle.t-particle_deer_collection[0].t0));
		draw(deer, scene);

		for(size_t k = 0; k<dementors.size(); k++){
			dementors[k].y0 += 0.27f;
			draw(dementors[k].visual, scene);
		}

	}
	else{
		arm.transform.rotate = vcl::rotation({1,0,0}, 0.0f);
		arm.transform.translate = {0.0f, -1.0f, 0.0f};

		for(size_t k = 0; k<dementors.size(); k++){
			if(dementors[k].y0>=2.0f)
				dementors[k].y0 -= 0.05f;
			draw(dementors[k].visual, scene);
		}


	}

	
	
	draw(forest.visual,scene);
	draw(lake, scene);
	draw_trees();
	//draw(deer, scene);

	//draw(mountain, scene);

	//character
	draw(head, scene);
	draw(body, scene);
	draw(arm, scene);
	draw(legs, scene);
	draw(hair1, scene);
	draw(hair2, scene);

	
	

	timer_particle.update();
	timer_spell.update();
	//if (timer_particle.event)
	//	particle_deer_collection.push_back(create_new_particle(timer_particle.t));

	for (size_t k = 0; k < particle_deer_collection.size(); ++k)
	{
		vec3 const p = compute_particle_position(particle_deer_collection[k], timer_particle.t);
		sphere.transform.translate = p;
		sphere.shading.color = particle_deer_collection[k].color;
		sphere.shading.alpha = 0.25 * ( 3 - (timer_particle.t - particle_deer_collection[k].t0 ) );
		sphere.transform.rotate = scene.camera.orientation();
		sphere.transform.scale = 0.3f;
		sphere.texture = opengl_texture_to_gpu(image_load_png("assets/sparcle.png"));
		draw(sphere, scene);
	}

	for (size_t k = 0; k < particle_spell_collection.size(); ++k)
	{
		vec3 const p = compute_particle_spell_position(particle_spell_collection[k], timer_spell.t);
		sphere2.transform.translate = p;
		//sphere.shading.color = particle_deer_collection[k].color;
		sphere2.shading.alpha = 0.25 * ( 3 - (timer_spell.t - particle_spell_collection[k].t0 ) );
		sphere2.transform.rotate = scene.camera.orientation();
		sphere.transform.scale = 0.5f;
		sphere2.shading.alpha = 0.5;
		draw(sphere2, scene);
	}


	for(size_t k = 0; k<dementors.size(); k++){
		float x0 = dementors[k].x0;
		float y0 = dementors[k].y0;
		float z0 = dementors[k].z0;
		float fi = dementors[k].fi;

		dementors[k].visual.transform.translate = {x0, 3.0f*cos(2*3.14f*timer.t + fi)+z0, y0};
		/*
		size_t const N_substeps = 5;
		for(size_t k_substep=0; k_substep<N_substeps; ++k_substep){
			compute_forces(dementors[k].forces, dementors[k].position, dementors[k].velocity, dementors[k].normal, 20.0f);
			numerical_integration(dementors[k].position, dementors[k].velocity, dementors[k].forces, m, dt);
			apply_constraints(dementors[k].position, dementors[k].velocity, dementors[k].positional_constraints, dementors[k].z0);	
			
		}

		dementors[k].visual.update_position(dementors[k].position.data);
		normal_per_vertex(dementors[k].position.data, dementors[k].triangle_connectivity, dementors[k].normal.data);
		dementors[k].visual.update_normal(dementors[k].normal.data);	
		*/
		draw(dementors[k].visual, scene);

	}



	glEnable(GL_BLEND);
	//glDepthMask(false);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	remove_old_particles(particle_deer_collection, timer_particle.t, 3.0f);
	remove_old_particles(particle_spell_collection, timer_spell.t, 4.0f);
}

void display_interface()
{
	ImGui::Checkbox("Display frame", &user.gui.display_frame);
	ImGui::Spacing(); ImGui::Spacing();

	ImGui::SliderFloat("Time", &timer.t, timer.t_min, timer.t_max, "%.2f s");
	ImGui::SliderFloat("Time Scale", &timer.scale, 0.05f, 2.0f, "%.2f s");


	
}


void window_size_callback(GLFWwindow* , int width, int height)
{
	glViewport(0, 0, width, height);
	float const aspect = width / static_cast<float>(height);
	scene.projection = projection_perspective(50.0f*pi/180.0f, aspect, 0.1f, 100.0f);
}


void mouse_move_callback(GLFWwindow* window, double xpos, double ypos)
{
	vec2 const  p1 = glfw_get_mouse_cursor(window, xpos, ypos);
	vec2 const& p0 = user.mouse_prev;

	glfw_state state = glfw_current_state(window);

	auto& camera = scene.camera;
	if(!user.cursor_on_gui){
		if(state.mouse_click_left && !state.key_ctrl)
			scene.camera.manipulator_rotate_trackball(p0, p1);
			// scene.camera.manipulator_rotate_trackball(p0 * vec2{ 1.0, 0.0 }, p1 * vec2{1.0, 0.0});
		if(state.mouse_click_left && state.key_ctrl)
			camera.manipulator_translate_in_plane(p1-p0);
			//camera.manipulator_translate_in_plane((p1 - p0) * (vec2{1.0, 0.0}));
		if(state.mouse_click_right)
			camera.manipulator_scale_distance_to_center( 0.1 * (p1-p0).y );
	}

	user.mouse_prev = p1;
}

void opengl_uniform(GLuint shader, scene_environment const& current_scene)
{
	opengl_uniform(shader, "projection", current_scene.projection);
	opengl_uniform(shader, "view", scene.camera.matrix_view());
	opengl_uniform(shader, "light", scene.light, false);
}



void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
	// TO DO
	if (key == GLFW_KEY_E){
		if (action == GLFW_PRESS || action == GLFW_REPEAT){
			levitate = true;

		}
		else{
			levitate = false;
		}
	}

}



