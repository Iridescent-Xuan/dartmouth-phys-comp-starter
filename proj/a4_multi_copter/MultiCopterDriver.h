#ifndef __MultiCopterDriver_h__
#define __MultiCopterDriver_h__
#include <random>
#include <iostream>
#include <fstream>
#include "Common.h"
#include "Driver.h"
#include "OpenGLMesh.h"
#include "OpenGLCommon.h"
#include "OpenGLWindow.h"
#include "OpenGLViewer.h"
#include "OpenGLParticles.h"
#include "Particles.h"
#include "MultiCopter.h"

template<int d> class MultiCopterDriver : public Driver, public OpenGLViewer
{
	using Base = Driver;
	double dt = .02;
	MultiCopter<d> copter;
	std::vector<Vector3> targets;
	std::vector<Vector3> waypoints;
	std::vector<Vector3> obstacles;

	OpenGLSegmentMesh* opengl_copter = nullptr;
	OpenGLSegmentMesh* opengl_circles = nullptr;
	std::vector<Vector3> circle_vtx;
	OpenGLParticles<Particles<3> >* opengl_targets = nullptr;
	OpenGLSegmentMesh* opengl_ground = nullptr;
	OpenGLSegmentMesh* opengl_obstacles = nullptr;

public:
	MultiCopterDriver() {

		//// default targets and obstacles
		targets = {
			Vector3(0.0, 0.0, -1.0),
			Vector3(1.0, 1.0, -1.2),
			Vector3(2.0, 0.0, -1.4),
			Vector3(1.0, -1.0, -1.2),
			Vector3(0.0, 0.0, -1.0),
			Vector3(0.0, 0.0, 0.0)
		};

		obstacles = {
			Vector3(1.,1.,0.),
			Vector3(2.,2.,0.),
			Vector3(3.,3.,0.)
		};

		//// -- TASK 3.5: Specify your own targets and obstacles here --
		// -- Your implementation starts --

		// -- Your implementation ends --


		// Create waypoints from targets.
		Vector3 last_target; last_target.setZero();
		const double max_dist_between_waypoints = 0.2;
		const int target_num = (int)targets.size();
		std::vector<Vector3> new_targets;
		for (const auto& next_target : targets) {
			const double target_dist = (next_target - last_target).norm();
			if (target_dist < max_dist_between_waypoints) {
				std::cout << "Warning: the target you set (" << next_target.transpose()
					<< ") is too close to its previous target (" << last_target.transpose()
					<< ") and will be ignored." << std::endl;
				continue;
			}

			const int waypoint_num = int(target_dist / max_dist_between_waypoints) + 1;
			const double waypoint_dist = target_dist / waypoint_num;
			for (int i = 0; i < waypoint_num; ++i)
				waypoints.push_back(last_target + (next_target - last_target) / target_dist * (i + 1) * waypoint_dist);

			last_target = next_target;
			new_targets.push_back(next_target);
		}

		new_targets.swap(targets);
	}

	virtual void Initialize(const int flag)
	{
		copter.Initialize(flag);

		const Vector3 dir = -Vector3::UnitZ();
		const double d1 = copter.arm_length / (double)(std::sqrt(2.0));
		copter.Add_Rotor(Vector3(d1, -d1, 0.0), dir);
		copter.Add_Rotor(Vector3(d1, d1, 0.0), dir);
		copter.Add_Rotor(Vector3(-d1, d1, 0.0), dir);
		copter.Add_Rotor(Vector3(-d1, -d1, 0.0), dir);

		////viewer initialization, initialize visualization data
		OpenGLViewer::Initialize();
		opengl_window->Set_Drone_View();
	}

	////synchronize simulation data to visualization data, called in OpenGLViewer::Initialize()
	virtual void Initialize_Data()
	{
		//// copter
		opengl_copter = Add_Interactive_Object<OpenGLSegmentMesh>();
		opengl_copter->mesh.elements.resize(2);
		opengl_copter->mesh.elements[0] = Vector2i(0, 2);
		opengl_copter->mesh.elements[1] = Vector2i(1, 3);
		opengl_copter->line_width = 6.f;
		*opengl_copter->mesh.vertices = copter.body_rotor_pos;

		opengl_copter->Set_Color(OpenGLColor(0.f, 1.f, 1.f, 1.f));
		opengl_copter->Set_Data_Refreshed();
		opengl_copter->Initialize();

		opengl_circles = Add_Interactive_Object<OpenGLSegmentMesh>();
		SegmentMesh<3>& mesh = opengl_circles->mesh;

		for (int i = 0; i < copter.body_rotor_pos.size(); i++) {
			Vector3 center = copter.body_rotor_pos[i];
			double r = .02f; int n = 16;
			double theta = 3.1415927f * 2.f / (double)n;
			int start = (*mesh.vertices).size();
			for (int j = 0; j < n; j++) {
				double angle = (double)j * theta;
				Vector3 p = center + Vector3(r * cos(angle), r * sin(angle), -0.005);
				circle_vtx.push_back(p);
				(*mesh.vertices).push_back(copter.World_Coord(p));
			}
			for (int j = 0; j < n - 1; j++) {
				mesh.elements.push_back(Vector2i(start + j, start + j + 1));
			}
			mesh.elements.push_back(Vector2i(start + n - 1, start));
		}

		opengl_circles->Set_Color(OpenGLColor(0.f, 1.f, 0.f, 1.f));

		opengl_circles->Set_Data_Refreshed();
		opengl_circles->Initialize();

		//// target points
		opengl_targets = Add_Interactive_Object<OpenGLParticles<Particles<3> > >();
		opengl_targets->particles.Resize(targets.size());
		for (int i = 0; i < targets.size(); i++) {
			opengl_targets->particles.X(i) = targets[i];
		}
		opengl_targets->Set_Data_Refreshed();
		opengl_targets->Initialize();

		//// ground
		opengl_ground = Add_Interactive_Object<OpenGLSegmentMesh>();
		SegmentMesh<3>& ground_mesh = opengl_ground->mesh;
		int n = 16;
		double dx = .5;
		Vector3 start = Vector3(1., 1., 0.) * (double)(-n / 2) * dx;
		for (int i = 0; i <= n; i++) {
			Vector3 p0 = start + Vector3::Unit(0) * (double)i * dx;
			Vector3 p1 = p0 + Vector3::Unit(1) * (double)n * dx;
			ground_mesh.Vertices().push_back(p0);
			ground_mesh.Vertices().push_back(p1);
			int s = ground_mesh.Vertices().size();
			ground_mesh.Elements().push_back(Vector2i(s - 2, s - 1));
		}
		for (int i = 0; i <= n; i++) {
			Vector3 p0 = start + Vector3::Unit(1) * (double)i * dx;
			Vector3 p1 = p0 + Vector3::Unit(0) * (double)n * dx;
			ground_mesh.Vertices().push_back(p0);
			ground_mesh.Vertices().push_back(p1);
			int s = ground_mesh.Vertices().size();
			ground_mesh.Elements().push_back(Vector2i(s - 2, s - 1));
		}

		opengl_ground->Set_Color(OpenGLColor(1.f, 1.f, 1.f, 1.f));
		opengl_ground->Set_Data_Refreshed();
		opengl_ground->Initialize();

		//// obstacles
		opengl_obstacles = Add_Interactive_Object<OpenGLSegmentMesh>();
		SegmentMesh<3>& obstacle_mesh = opengl_obstacles->mesh;
		for (int i = 0; i < obstacles.size(); i++) {
			obstacle_mesh.Vertices().push_back(obstacles[i]);
			obstacle_mesh.Vertices().push_back(obstacles[i] + Vector3(0., 0., -1.));
			int s = obstacle_mesh.Vertices().size();
			obstacle_mesh.Elements().push_back(Vector2i(s - 2, s - 1));
		}
		opengl_obstacles->line_width = 6.f;
		opengl_obstacles->Set_Color(OpenGLColor(1.f, 1.f, 0.f, 1.f));
		opengl_obstacles->Set_Data_Refreshed();
		opengl_obstacles->Initialize();
	}

	void Sync_Simulation_And_Visualization_Data()
	{
		for (int i = 0; i < copter.body_rotor_pos.size(); i++) {
			(*opengl_copter->mesh.vertices)[i] = copter.World_Coord(copter.body_rotor_pos[i]);
		}
		opengl_copter->Set_Data_Refreshed();

		for (int i = 0; i < (*opengl_circles->mesh.vertices).size(); i++) {
			(*opengl_circles->mesh.vertices)[i] = copter.World_Coord(circle_vtx[i]);
		}
		opengl_circles->Set_Data_Refreshed();

	}

	////update simulation and visualization for each time step
	virtual void Toggle_Next_Frame()
	{
		static double t = 0.0;
		// We give each waypoint 1 second.
		int idx = (int)(t / 1.0);
		if (idx > (int)waypoints.size() - 1) idx = (int)waypoints.size() - 1;
		copter.Advance(dt, waypoints[idx]);
		Sync_Simulation_And_Visualization_Data();
		OpenGLViewer::Toggle_Next_Frame();
		t += dt;
	}

	virtual void Run()
	{
		OpenGLViewer::Run();
	}
};
#endif