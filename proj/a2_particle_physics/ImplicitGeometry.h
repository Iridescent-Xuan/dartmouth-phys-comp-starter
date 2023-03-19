#ifndef __ImplicitGeometry_h__
#define __ImplicitGeometry_h__

template<int d> class ImplicitGeometry
{
	using VectorD = Vector<double, d>;
public:
	virtual double Phi(const VectorD& pos) const { return 0.; }
	virtual VectorD Normal(const VectorD& pos) const { return VectorD::Zero(); }
};

template<int d> class Bowl : public ImplicitGeometry<d>
{
	using VectorD = Vector<double, d>;
public:
	VectorD center;
	double radius;
	Bowl(VectorD _center = VectorD::Zero(), double _radius = 1.) :center(_center), radius(_radius) {}
	virtual double Phi(const VectorD& pos) const { return radius - (pos - center).norm(); }
	virtual VectorD Normal(const VectorD& pos) const { return (center - pos).normalized(); }
};

template<int d> class Sphere : public ImplicitGeometry<d>
{
	using VectorD = Vector<double, d>;
public:
	VectorD center;
	double radius;
	Sphere(VectorD _center = VectorD::Zero(), double _radius = 1.) :center(_center), radius(_radius) {}
	virtual double Phi(const VectorD& pos) const { return (pos - center).norm() - radius; }
	virtual VectorD Normal(const VectorD& pos) const { return (pos - center).normalized(); }
};

#endif
