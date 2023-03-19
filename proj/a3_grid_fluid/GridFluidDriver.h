#ifndef __GridFluidDriver_h__
#define __GridFluidDriver_h__
#include <random>
#include "Common.h"
#include "Driver.h"
#include "Particles.h"
#include "OpenGLMesh.h"
#include "OpenGLCommon.h"
#include "OpenGLWindow.h"
#include "OpenGLViewer.h"
#include "OpenGLMarkerObjects.h"
#include "OpenGLParticles.h"

class GridFluidDriver : public Driver, public OpenGLViewer
{
	using Base = Driver;
	double dt = .02;
	GridFluid fluid;
	double v_rescale = (double).05;
	bool add_particle = false;

	OpenGLSegmentMesh* opengl_vectors = nullptr;							////vector field for velocity
	OpenGLPolygon* opengl_polygon = nullptr;								////a rectangle for domain boundary
	std::vector<OpenGLSolidCircle*> opengl_circles;							////passive particles
	Hashset<int> invis_particles;
	OpenGLColoredTriangleMesh* opengl_mesh = nullptr;						////density field

	//////////////////////////////////////////////////////////////////////////
	////specify the following flags to control visualization
	bool draw_velocity = true;
	bool draw_density = true;
	bool draw_particles = true;
	//////////////////////////////////////////////////////////////////////////

public:
	virtual void Initialize()
	{
		fluid.Initialize();
		fluid.Initialize_Visualization_Particles();
		OpenGLViewer::Initialize();
	}

	////initialize visualization data, called in OpenGLViewer::Initialize()
	virtual void Initialize_Data()
	{
		////initialize vector field (write all vectors as segments)
		opengl_vectors = Add_Interactive_Object<OpenGLSegmentMesh>();
		opengl_vectors->mesh.elements.resize(fluid.node_num);
		opengl_vectors->mesh.vertices->resize(fluid.node_num * 2);

		for (int i = 0; i < fluid.particles.Size(); i++) {
			Add_Solid_Circle(i);
		}

		for (int i = 0; i < fluid.node_num; i++) {
			opengl_vectors->mesh.elements[i] = Vector2i(i * 2, i * 2 + 1);
			(*opengl_vectors->mesh.vertices)[i * 2] = V3(fluid.grid.Node(fluid.grid.Node_Coord(i)));
			Vector2 pos2 = fluid.grid.Node(fluid.grid.Node_Coord(i)) + fluid.u[i] * v_rescale;
			(*opengl_vectors->mesh.vertices)[i * 2 + 1] = V3(pos2);
		}
		opengl_vectors->Set_Data_Refreshed();
		opengl_vectors->Initialize();

		////initialize polygon
		opengl_polygon = Add_Interactive_Object<OpenGLPolygon>();
		opengl_polygon->vtx.push_back(Vector3::Zero());
		opengl_polygon->vtx.push_back(Vector3::Unit(0) * (double)2);
		opengl_polygon->vtx.push_back(Vector3::Unit(0) * (double)2 + Vector3::Unit(1) * (double)1);
		opengl_polygon->vtx.push_back(Vector3::Unit(1) * (double)1);
		Set_Color(opengl_polygon, OpenGLColor(.0, 1., 1., 1.));
		Set_Line_Width(opengl_polygon, 4.f);
		opengl_polygon->Set_Data_Refreshed();
		opengl_polygon->Initialize();

		opengl_mesh = Add_Interactive_Object<OpenGLColoredTriangleMesh>();

		int cell_num = fluid.grid.cell_counts.prod();
		for (int i = 0; i < cell_num; i++) {
			Vector2i cell = fluid.grid.Cell_Coord(i);
			Vector2 pos1 = fluid.grid.Node(cell);
			Vector2 pos2 = pos1 + Vector2::Unit(0) * fluid.grid.dx;
			Vector2 pos3 = pos1 + Vector2::Unit(0) * fluid.grid.dx + Vector2::Unit(1) * fluid.grid.dx;
			Vector2 pos4 = pos1 + Vector2::Unit(1) * fluid.grid.dx;

			int n = (int)opengl_mesh->mesh.Vertices().size();
			opengl_mesh->mesh.Vertices().push_back(V3(pos1));
			opengl_mesh->mesh.Vertices().push_back(V3(pos2));
			opengl_mesh->mesh.Vertices().push_back(V3(pos3));
			opengl_mesh->colors.push_back(0.);
			opengl_mesh->colors.push_back(0.);
			opengl_mesh->colors.push_back(0.);
			opengl_mesh->mesh.Elements().push_back(Vector3i(n, n + 1, n + 2));

			n = (int)opengl_mesh->mesh.Vertices().size();
			opengl_mesh->mesh.Vertices().push_back(V3(pos1));
			opengl_mesh->mesh.Vertices().push_back(V3(pos3));
			opengl_mesh->mesh.Vertices().push_back(V3(pos4));
			opengl_mesh->colors.push_back(0.);
			opengl_mesh->colors.push_back(0.);
			opengl_mesh->colors.push_back(0.);
			opengl_mesh->mesh.Elements().push_back(Vector3i(n, n + 1, n + 2));
		}

		opengl_mesh->Set_Data_Refreshed();
		opengl_mesh->Initialize();
	}

	////synchronize simulation data to visualization data
	void Sync_Simulation_And_Visualization_Data()
	{
		////velocity visualization
		if (draw_velocity) {
			for (int i = 0; i < fluid.node_num; i++) {
				Vector2 pos2 = fluid.grid.Node(fluid.grid.Node_Coord(i)) + fluid.u[i] * v_rescale;
				(*opengl_vectors->mesh.vertices)[i * 2 + 1] = V3(pos2);
			}
			opengl_vectors->Set_Data_Refreshed();
		}

		////particle visualization
		if (draw_particles) {
			fluid.Update_Visualization_Particles(dt);
			for (int i = 0; i < fluid.particles.Size(); i++) {
				bool outside = false;
				double epsilon = fluid.grid.dx * (double).5;
				for (int j = 0; j < 2; j++) {
					if (fluid.particles.X(i)[j]<fluid.grid.domain_min[j] + epsilon ||
						fluid.particles.X(i)[j]>fluid.grid.domain_max[j] - epsilon) {
						outside = true; break;
					}
				}
				if (outside) {
					opengl_circles[i]->visible = false;
					opengl_circles[i]->Set_Data_Refreshed();
					invis_particles.insert(i);
					fluid.particles.I(i) = -1;
					continue;
				}

				auto opengl_circle = opengl_circles[i];
				opengl_circle->pos = V3(fluid.particles.X(i));
				opengl_circle->pos[2] = (double).05;
				opengl_circle->Set_Data_Refreshed();
			}
		}

		////density visualization
		if (draw_density) {
			int cell_num = fluid.grid.cell_counts.prod();
			for (int i = 0; i < cell_num; i++) {
				int idx = i * 6;
				Vector2i cell = fluid.grid.Cell_Coord(i);
				Vector2 pos1 = fluid.grid.Node(cell);
				double den1 = fluid.Interpolate(fluid.smoke_den, pos1);

				Vector2 pos2 = pos1 + Vector2::Unit(0) * fluid.grid.dx;
				double den2 = fluid.Interpolate(fluid.smoke_den, pos2);

				Vector2 pos3 = pos1 + Vector2::Unit(0) * fluid.grid.dx + Vector2::Unit(1) * fluid.grid.dx;
				double den3 = fluid.Interpolate(fluid.smoke_den, pos3);

				Vector2 pos4 = pos1 + Vector2::Unit(1) * fluid.grid.dx;
				double den4 = fluid.Interpolate(fluid.smoke_den, pos4);

				opengl_mesh->colors[idx] = den1;
				opengl_mesh->colors[idx + 1] = den2;
				opengl_mesh->colors[idx + 2] = den3;

				opengl_mesh->colors[idx + 3] = den1;
				opengl_mesh->colors[idx + 4] = den3;
				opengl_mesh->colors[idx + 5] = den4;
			}
			opengl_mesh->Set_Data_Refreshed();
		}
	}

	////update simulation and visualization for each time step
	virtual void Toggle_Next_Frame()
	{
		fluid.Advance(dt);
		Sync_Simulation_And_Visualization_Data();
		OpenGLViewer::Toggle_Next_Frame();
	}

	virtual void Run()
	{
		OpenGLViewer::Run();
	}

	////User interaction
	virtual bool Mouse_Drag(int x, int y, int w, int h)
	{
		if (!add_particle)return false;
		Vector3f win_pos = opengl_window->Project(Vector3f::Zero());
		Vector3f pos = opengl_window->Unproject(Vector3f((float)x, (float)y, win_pos[2]));
		Vector2 p_pos; for (int i = 0; i < 2; i++)p_pos[i] = (double)pos[i];
		fluid.src_pos = p_pos;
		Add_Source_Particle(p_pos);
		return true;
	}

	virtual bool Mouse_Click(int left, int right, int mid, int x, int y, int w, int h)
	{
		if (left != 1 && left != -1) { return false; }
		if (left == -1 && add_particle) {
			add_particle = false;
			fluid.src_pos = Vector2::Ones() * -1;	////turn off source
			return true;
		}

		Vector3f win_pos = opengl_window->Project(Vector3f::Zero());
		Vector3f pos = opengl_window->Unproject(Vector3f((float)x, (float)y, win_pos[2]));
		Vector2 p_pos; for (int i = 0; i < 2; i++)p_pos[i] = (double)pos[i];
		fluid.src_pos = p_pos;
		Add_Source_Particle(p_pos);
		add_particle = true;
		return true;
	}

	////Keyboard interaction
	virtual void Initialize_Common_Callback_Keys()
	{
		OpenGLViewer::Initialize_Common_Callback_Keys();
		Bind_Callback_Key('v', &Keyboard_Event_V_Func, "press v");
		Bind_Callback_Key('2', &Keyboard_Event_V_Func, "press 2");
	}

	virtual void Keyboard_Event_V()
	{
		draw_velocity = !draw_velocity;
		opengl_vectors->visible = !opengl_vectors->visible;
		opengl_vectors->Set_Data_Refreshed();
	}
	Define_Function_Object(GridFluidDriver, Keyboard_Event_V);

	virtual void Keyboard_Event_D()
	{
		draw_density = !draw_density;
		opengl_mesh->visible = !opengl_mesh->visible;
		opengl_mesh->Set_Data_Refreshed();
	}
	Define_Function_Object(GridFluidDriver, Keyboard_Event_D);

	//////////////////////////////////////////////////////////////////////////
	////User interaction for manipulating sources
	void Add_Source_Particle(Vector2 p_pos)
	{
		double rx = .1 * static_cast<float>(rand() % 1000) / 1000. - .05;
		double ry = .1 * static_cast<float>(rand() % 1000) / 1000. - .05;
		p_pos[0] += rx; p_pos[1] += ry;
		if (!invis_particles.empty()) {
			int p = (*invis_particles.begin());
			fluid.particles.X(p) = p_pos;
			fluid.particles.I(p) = 0;
			opengl_circles[p]->pos = V3(fluid.particles.X(p));
			opengl_circles[p]->Update_Model_Matrix();
			opengl_circles[p]->Set_Data_Refreshed();
			opengl_circles[p]->visible = true;
			invis_particles.erase(p);
		}
		else {
			Add_Particle(p_pos);
			Add_Solid_Circle(fluid.particles.Size() - 1);
		}
	}

	void Add_Particle(Vector2 pos)
	{
		int i = fluid.particles.Add_Element();	////return the last element's index
		fluid.particles.X(i) = pos;
		fluid.particles.C(i) = (double)(rand() % 2000 - 1000) / (double)1000;	////particle vorticity, a random number between [-1,1]
	}

	void Add_Solid_Circle(const int i)
	{
		OpenGLColor c(1.f, 0.68f, 0.26f, 1.f);
		auto opengl_circle = Add_Interactive_Object<OpenGLSolidCircle>();
		opengl_circle->visible = false;
		opengl_circles.push_back(opengl_circle);
		opengl_circle->pos = V3(fluid.particles.X(i));
		opengl_circle->radius = (double).01;
		opengl_circle->color = c;
		opengl_circle->Initialize();
		opengl_circle->Update_Model_Matrix();
		opengl_circle->Set_Data_Refreshed();
		opengl_circle->visible = true;
	}

protected:
	////Helper function to convert a vector to 3d, for c++ template
	Vector3 V3(const Vector2& v2) { return Vector3(v2[0], v2[1], .0); }
	Vector3 V3(const Vector3& v3) { return v3; }
};
#endif