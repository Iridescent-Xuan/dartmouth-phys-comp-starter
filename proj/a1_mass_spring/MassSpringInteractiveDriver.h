//#####################################################################
// Mass spring driver
// Dartmouth COSC 89.18/189.02: Computational Methods for Physical Systems, Assignment starter code
// Contact: Bo Zhu (bo.zhu@dartmouth.edu)
//#####################################################################
#ifndef __MassSpringInteractiveDriver_h__
#define __MassSpringInteractiveDriver_h__
#include <memory>
#include "Common.h"
#include "Mesh.h"
#include "InClassDemoDriver.h"
#include "SoftBodyMassSpring.h"

class MassSpringInteractivDriver : public InClassDemoDriver
{
	using Base = Driver;
public:
	SoftBodyMassSpring soft_body;
	const double dt = (double).02;

	////visualization data
	Segments segments;
	std::vector<Point> points;

	////initialize simulation data and its visualizations
	virtual void Initialize_Data()
	{
		Initialize_Simulation_Data();

		segments.Initialize(this);
		segments.Sync_Data(soft_body.particles.XRef(), soft_body.springs);

		int n = soft_body.particles.Size();
		points.resize(n);
		for (int i = 0; i < n; i++) {
			points[i].Set_Radius(.02);
			points[i].Initialize(this);
			points[i].Sync_Data(soft_body.particles.X(i));
		}
	}

	////advance simulation timesteps
	virtual void Advance(const double dt)
	{
		soft_body.Advance(dt);
	}

	////update simulation data to its visualization counterparts
	virtual void Sync_Simulation_And_Visualization_Data()
	{
		segments.Sync_Data(soft_body.particles.XRef());
		int n = soft_body.particles.Size();
		for (int i = 0; i < n; i++) {
			points[i].Sync_Data(soft_body.particles.X(i));
		}
	}

	////here we initialize three tests for rod, cloth, and beam
	virtual void Initialize_Simulation_Data()
	{
		switch (test) {
		case 1: {	////1d rod
			////initialize spring vertices
			double length = (double)1; int n = 8; double dx = length / (double)n;
			soft_body.particles.Resize(n);
			for (int i = 0; i < n; i++) {
				soft_body.particles.X(i) = Vector3::Unit(0) * (double)i * dx;
				soft_body.particles.M(i) = (double)1;
			}
			////initialize springs
			for (int i = 0; i < n - 1; i++) {
				Vector2i s(i, i + 1);
				soft_body.springs.push_back(s);
			}
			////set boundary conditions
			soft_body.Set_Boundary_Node(0);
		}break;
		case 2: {	////2d cloth
			////create a cloth mesh
			double length = (double)1; int width = 4 * scale; int height = 6 * scale; double step = length / (double)width;
			TriangleMesh<3> cloth_mesh;
			Build_Cloth_Mesh(width, height, step, &cloth_mesh, 0, 2);
			int n = (int)cloth_mesh.Vertices().size();
			std::vector<Vector2i> edges; Get_Mesh_Edges(cloth_mesh, edges);

			////copy cloth mesh vertices to spring particles 
			soft_body.particles.Resize(n);
			for (int i = 0; i < n; i++) {
				soft_body.particles.X(i) = cloth_mesh.Vertices()[i];
				soft_body.particles.M(i) = (double)1;
			}
			////copy cloth mesh edges to springs
			soft_body.springs = edges;

			////set boundary conditions
			soft_body.Set_Boundary_Node(0);
			soft_body.Set_Boundary_Node(width - 1);
		}break;
		case 3: {	////3d volumetric beam
			int n = 4 * scale; double dx = (double)1 / (double)n;
			Build_Beam_Particles_And_Springs(soft_body.particles, soft_body.springs, n, dx);
			for (int i = 0; i < 4; i++)soft_body.Set_Boundary_Node(i);
		}break;

			//////////////////////////////////////////////////////////////////////////
			////YOUR IMPLEMENTATION (TASK 2: OPTION 1): simulate a single hair strand
		case 4: {
			soft_body.Initialize_Hair_Strand();
		}break;
		}

		soft_body.Initialize();
	}

	//////////////////////////////////////////////////////////////////////////
	////These helper functions are all for creating meshes
protected:
	////Helper functions
	void Build_Cloth_Mesh(const int cell_num_0, const int cell_num_1, const double dx, TriangleMesh<3>* mesh, int axis_0 = 0, int axis_1 = 1)
	{
		mesh->elements.resize(2 * (cell_num_0 - 1) * (cell_num_1 - 1)); int t = 0;
		for (int i = 1; i <= cell_num_0 - 1; i++)for (int j = 1; j <= cell_num_1 - 1; j++) { // counterclockwise node ordering
			if (i % 2) { mesh->elements[t++] = Vector3i(i + cell_num_0 * (j - 1), i + 1 + cell_num_0 * (j - 1), i + cell_num_0 * j); mesh->elements[t++] = Vector3i(i + 1 + cell_num_0 * (j - 1), i + 1 + cell_num_0 * j, i + cell_num_0 * j); }
			else { mesh->elements[t++] = Vector3i(i + cell_num_0 * (j - 1), i + 1 + cell_num_0 * (j - 1), i + 1 + cell_num_0 * j); mesh->elements[t++] = Vector3i(i + cell_num_0 * (j - 1), i + 1 + cell_num_0 * j, i + cell_num_0 * j); }
		}
		for (size_type i = 0; i < mesh->elements.size(); i++) {
			mesh->elements[i] -= Vector3i::Ones();
			/*swap y and z*/int tmp = mesh->elements[i][1]; mesh->elements[i][1] = mesh->elements[i][2]; mesh->elements[i][2] = tmp;
		}
		for (int j = 0; j < cell_num_1; j++)for (int i = 0; i < cell_num_0; i++) { Vector3 pos = Vector3::Zero(); pos[axis_0] = (double)i * dx; pos[axis_1] = (double)j * dx; mesh->Vertices().push_back(pos); }
	}

	void Get_Mesh_Edges(const TriangleMesh<3>& mesh, std::vector<Vector2i>& edges)
	{
		Hashset<Vector2i> edge_hashset; ArrayF<Vector2i, 6> element_edges;
		for (const auto& vtx : mesh.elements) {
			edge_hashset.insert(Sorted(Vector2i(vtx[0], vtx[1])));
			edge_hashset.insert(Sorted(Vector2i(vtx[1], vtx[2])));
			edge_hashset.insert(Sorted(Vector2i(vtx[2], vtx[0])));
		}
		for (const auto& edge : edge_hashset)edges.push_back(edge);
	}

	void Build_Beam_Particles_And_Springs(Particles<3>& particles, std::vector<Vector2i>& edges, int n, double dx, Vector3 pos = Vector3::Zero())
	{
		particles.Resize(n * 4);
		for (int i = 0; i < particles.Size(); i++) {
			particles.M(i) = (double)1;
		}
		for (int i = 0; i < n; i++) {
			particles.X(i * 4) = pos + Vector3(dx * (double)i, (double)0, (double)0);
			particles.X(i * 4 + 1) = pos + Vector3(dx * (double)i, (double)0, (double)dx);
			particles.X(i * 4 + 2) = pos + Vector3(dx * (double)i, (double)dx, (double)0);
			particles.X(i * 4 + 3) = pos + Vector3(dx * (double)i, (double)dx, (double)dx);
			edges.push_back(Vector2i(i * 4, i * 4 + 1));
			edges.push_back(Vector2i(i * 4 + 1, i * 4 + 3));
			edges.push_back(Vector2i(i * 4 + 3, i * 4 + 2));
			edges.push_back(Vector2i(i * 4 + 2, i * 4));
			if (i < n - 1) {
				edges.push_back(Vector2i(i * 4, i * 4 + 4));
				edges.push_back(Vector2i(i * 4 + 1, i * 4 + 5));
				edges.push_back(Vector2i(i * 4 + 2, i * 4 + 6));
				edges.push_back(Vector2i(i * 4 + 3, i * 4 + 7));

				edges.push_back(Vector2i(i * 4, i * 4 + 7));
				edges.push_back(Vector2i(i * 4 + 1, i * 4 + 6));
				edges.push_back(Vector2i(i * 4 + 2, i * 4 + 5));
				edges.push_back(Vector2i(i * 4 + 3, i * 4 + 4));
			}
		}
	}

	Vector2i Sorted(const Vector2i& v) { return v[0] > v[1] ? v : Vector2i(v[1], v[0]); }
};
#endif