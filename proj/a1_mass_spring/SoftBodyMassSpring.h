//////////////////////////////////////////////////////////////////////////
//// Dartmouth Physical Computing Programming Assignment 1: Mass Spring
//// Author: TODO: PUT YOUR NAME HERE
////////////////////////////////////////////////////////////////////////// 

#ifndef __SoftBodyMassSpring_h__
#define __SoftBodyMassSpring_h__
#include "Common.h"
#include "Particles.h"

class SoftBodyMassSpring
{
public:
	////Spring parameters
	Particles<3> particles;								//// The particle system. Each particle has the attributes of position, velocity, mass, and force. Read Particles.h in src to know the details.
	std::vector<Vector2i> springs;						//// Each element in springs stores a Vector2i for the indices of the two end points.
	std::vector<double> rest_length;					//// Each element in rest_length specifies the rest length of the spring
	std::vector<double> ks;								//// Each element in ks specifies the spring stiffness of a spring
	std::vector<double> kd;								//// Each element in kd specifies the damping coefficient of a spring

	////Boundary nodes
	std::unordered_map<int, Vector3> boundary_nodes;		//// boundary_notes stores the mapping from node index to its specified velocity. E.g., a fixed node will have a zero velocity.

	////Body force
	Vector3 g = Vector3::Unit(1) * (double)-1.;			//// gravity

	enum class TimeIntegration { ExplicitEuler, ImplicitEuler } time_integration = TimeIntegration::ImplicitEuler;	//// set to ExplicitEuler by default; change it to ImplicitEuler when you work on Task 2 (Option 2)

	////Implicit time integration
	SparseMatrixT K;
	VectorX u, b;

	virtual void Initialize()
	{
		////Initialize default spring parameters for standard tests
		double ks_0 = (double)1, kd_0 = (double)1;
		switch (time_integration) {
		case TimeIntegration::ExplicitEuler: {
			ks_0 = (double)5e2;
			kd_0 = (double)1e1;
		}break;
		case TimeIntegration::ImplicitEuler: {
			ks_0 = (double)1e5;
			kd_0 = (double)1e1;
		}break;
		}

		////Allocate arrays for springs and parameters
		rest_length.resize(springs.size());
		for (int i = 0; i < (int)springs.size(); i++) {
			const Vector2i& s = springs[i];
			rest_length[i] = (particles.X(s[0]) - particles.X(s[1])).norm();
		}
		ks.resize(springs.size(), ks_0);
		kd.resize(springs.size(), kd_0);

		////Allocate sparse matrix if using implicit time integration 
		////This function needs to be called for only once since the mesh doesn't change during the simulation)
		if (time_integration == TimeIntegration::ImplicitEuler)
			Initialize_Implicit_K_And_b();
	}

	virtual void Advance(const double dt)
	{
		switch (time_integration) {
		case TimeIntegration::ExplicitEuler:
			Advance_Explicit_Euler(dt); break;
		case TimeIntegration::ImplicitEuler:
			Advance_Implicit_Euler(dt); break;
		}
	}

	////Set boundary nodes
	void Set_Boundary_Node(const int p, const Vector3 v = Vector3::Zero()) { boundary_nodes[p] = v; }

	bool Is_Boundary_Node(const int p) { return boundary_nodes.find(p) != boundary_nodes.end(); }

	//////////////////////////////////////////////////////////////////////////
	//// P1 TASK: explicit Euler integration and spring force calculation

	//////////////////////////////////////////////////////////////////////////
	//// YOUR IMPLEMENTATION (TASK 1): explicit Euler time integration 
	void Advance_Explicit_Euler(const double dt)
	{
		//// Step 0: Clear the force on each particle (already done for you)
		Clear_Force();

		//// Step 1: add a body force to each particle
		Apply_Body_Force(dt);

		//// Step 2: calculate the spring force for each spring and add it to the two connecting particles (Hint: you may want to implement the function Spring_Force and call it in your for-loop)
		Apply_Spring_Force(dt);

		//// Step 3: enforce the boundary conditions (traversing the particles in the unordered_set, set its velocity to be the value read from the unordered_set, and its force to be zero)
		Enforce_Boundary_Condition();

		//// Step 4: time integration by updating the particle velocities according to their forces and updating particle positions according to the positions
		Time_Integration(dt);
	}

	void Clear_Force()
	{
		for (int i = 0; i < particles.Size(); i++) { particles.F(i) = Vector3::Zero(); }
	}

	void Apply_Body_Force(const double dt)
	{
		/* Your implementation start */

		// TODO: what is @dt used for?
		for (int i = 0; i < particles.Size(); i++) { particles.F(i) += particles.M(i) * g; }

		/* Your implementation end */
	}

	void Apply_Spring_Force(const double dt)
	{
		/* Your implementation start */

		assert(springs.size() == kd.size());
		assert(springs.size() == ks.size());

		for (int i = 0; i < springs.size(); ++i) {
			int pi = springs[i][0];
			int pj = springs[i][1];

			particles.F(pi) += Spring_Force(i);
			particles.F(pj) -= Spring_Force(i);
		}

		/* Your implementation end */
	}

	void Enforce_Boundary_Condition()
	{
		/* Your implementation start */

		for (auto i : boundary_nodes) {
			particles.V(i.first) = i.second;
			particles.F(i.first) = Vector3::Zero();
		}

		/* Your implementation end */
	}

	void Time_Integration(const double dt)
	{
		/* Your implementation start */

		for (int i = 0; i < particles.Size(); ++i) {
			particles.V(i) += dt * particles.F(i) / particles.M(i);
			particles.X(i) += dt * particles.V(i);
		}

		/* Your implementation end */
	}

	Vector3 Spring_Force(const int spring_index)
	{
		//// This is an auxiliary function to compute the spring force f=f_s+f_d for the spring with spring_index. 
		//// You may want to call this function in Apply_Spring_Force

		/* Your implementation start */

		Vector2i spring = springs[spring_index];
		Vector3d l_ij = particles.X(spring[1]) - particles.X(spring[0]);
		Vector3 f_s = ks[spring_index] * (l_ij.norm() - rest_length[spring_index]) * l_ij.normalized();
		Vector3 f_d = kd[spring_index] * (particles.V(spring[1]) - particles.V(spring[0])).dot(l_ij.normalized()) * l_ij.normalized();

		/* Your implementation end */

		return f_s + f_d;	////REPLACE this line with your own implementation
	}

	//////////////////////////////////////////////////////////////////////////
	//// TASK 2 (OPTION 1): creating bending springs for hair strand simulation
	void Initialize_Hair_Strand()
	{
		//// You need to initialize a hair model by setting the springs and particles to simulate human hair.
		//// A key component for a hair simulator is the bending spring (e.g., by connecting particles with a certain index offset).
		//// Think about how to realize these bending and curly effects with the explicit spring model you have implemented in TASK 1.
		//// You may also want to take a look at the function Initialize_Simulation_Data() in MassSpringInteractiveDriver.h for the model initialization.

		/* Your implementation start */

		/* Your implementation end */
	}

	//////////////////////////////////////////////////////////////////////////
	//// TASK 2 (OPTION 2): implicit time integration for inextensible cloth simulation
	//// The rest part of this file is all for this task.

	////Construct K, step 1: initialize the matrix structure 
	void Initialize_Implicit_K_And_b()
	{
		int n = 3 * particles.Size();
		K.resize(n, n); u.resize(n); u.fill((double)0); b.resize(n); b.fill((double)0);
		std::vector<TripletT> elements;
		for (int s = 0; s < (int)springs.size(); s++) {
			int i = springs[s][0]; int j = springs[s][1];
			Add_Block_Triplet_Helper(i, i, elements);
			Add_Block_Triplet_Helper(i, j, elements);
			Add_Block_Triplet_Helper(j, i, elements);
			Add_Block_Triplet_Helper(j, j, elements);
		}
		K.setFromTriplets(elements.begin(), elements.end());
		K.makeCompressed();
	}

	////Construct K, step 2: fill nonzero elements in K
	void Update_Implicit_K_And_b(const double dt)
	{
		////Clear K and b
		K.setZero();
		b.fill((double)0);

		/* Your implementation start */

		for (int i = 0; i < particles.Size(); ++i) {
			// for K
			// M
			for (int d = 0; d < 3; ++d) {
				K.coeffRef(i * 3 + d, i * 3 + d) = particles.M(i);
			}

			// for b
			// M * v^n + dt * f^n
			Set_Block(b, i, particles.M(i) * particles.V(i) + dt * particles.F(i));
		}

		for (int i = 0; i < springs.size(); ++i) {
			int pi = springs[i][0];
			int pj = springs[i][1];

			Matrix3 Ks, Kd;
			Compute_Ks_Block(i, Ks);
			Compute_Kd_Block(i, Kd);

			// for K
			// -dt * D
			Add_Block_Helper(K, pi, pj, -dt * Kd);
			// -dt^2 * K
			Add_Block_Helper(K, pi, pj, -dt * dt * Ks);

			// for d
			Add_Block(b, pi, -dt * Kd * particles.V(pi));
			Add_Block(b, pj, dt * Kd * particles.V(pj));
		}

		/* Your implementation end */
	}

	////Construct K, step 2.1: compute spring force derivative
	void Compute_Ks_Block(const int s, Matrix3& Ks)
	{
		/* Your implementation start */

		int pi = springs[s][0];
		int pj = springs[s][1];
		Vector3d x_ji = particles.X(pj) - particles.X(pi);

		Ks = (rest_length[s] / x_ji.norm() - 1) * Matrix3d::Identity();
		Ks -= rest_length[s] * x_ji * x_ji.transpose() / pow(x_ji.norm(), 3);
		Ks *= ks[s];

		/* Your implementation end */
	}

	////Construct K, step 2.2: compute damping force derivative
	void Compute_Kd_Block(const int s, Matrix3& Kd)
	{
		/* Your implementation start */

		int pi = springs[s][0];
		int pj = springs[s][1];
		Vector3d x_ji = particles.X(pj) - particles.X(pi);
		Kd = -1.0 * x_ji * x_ji.transpose() / pow(x_ji.norm(), 2);
		Kd *= kd[s];

		/* Your implementation end */
	}

	////Implicit Euler time integration
	void Advance_Implicit_Euler(const double dt)
	{
		//// clear force
		Clear_Force();
		//// add a body force to each particle as in explicit Euler
		Apply_Body_Force(dt);
		//// calculate the spring force as in explicit Euler
		Apply_Spring_Force(dt);
		//// enforce boundary condition
		Enforce_Boundary_Condition();

		Update_Implicit_K_And_b(dt);

		for (int i = 0; i < particles.Size(); i++) {
			for (int j = 0; j < 3; j++)u[i * 3 + j] = particles.V(i)[j];
		}	////set initial guess to be the velocity from the last time step

		SparseSolver::CG(K, u, b);	////solve Ku=b using Conjugate Gradient

		for (int i = 0; i < particles.Size(); i++) {
			Vector3 v; for (int j = 0; j < 3; j++)v[j] = u[i * 3 + j];
			particles.V(i) = v;
			particles.X(i) += particles.V(i) * dt;
		}
	}

	////Hint: you may want to use these functions when assembling your implicit matrix
	////Add block nonzeros to sparse matrix elements (for initialization)
	void Add_Block_Triplet_Helper(const int i, const int j, std::vector<TripletT>& elements)
	{
		for (int ii = 0; ii < 3; ii++)for (int jj = 0; jj < 3; jj++)elements.push_back(TripletT(i * 3 + ii, j * 3 + jj, (double)0));
	}

	////Add block Ks to K_ij
	void Add_Block_Helper(SparseMatrixT& K, const int i, const int j, const Matrix3& Ks)
	{
		SparseFunc::Add_Block<3, Matrix3>(K, i, i, Ks);
		SparseFunc::Add_Block<3, Matrix3>(K, j, j, Ks);
		if (!Is_Boundary_Node(i) && !Is_Boundary_Node(j)) {
			SparseFunc::Add_Block<3, Matrix3>(K, i, j, -Ks);
			SparseFunc::Add_Block<3, Matrix3>(K, j, i, -Ks);
		}
	}

	////Set block values on a vector
	void Set_Block(VectorX& b, const int i, const Vector3& bi)
	{
		for (int ii = 0; ii < 3; ii++)b[i * 3 + ii] = bi[ii];
	}

	////Add block values to a vector
	void Add_Block(VectorX& b, const int i, const Vector3& bi)
	{
		for (int ii = 0; ii < 3; ii++)b[i * 3 + ii] += bi[ii];
	}
};

#endif