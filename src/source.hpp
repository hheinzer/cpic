#ifndef SOURCE_HPP
#define SOURCE_HPP

#include <iostream>
#include <Eigen/Dense>
#include "domain.hpp"
#include "species.hpp"

class Source {
	public:
		using Vector3d = Eigen::Vector3d;

		virtual ~Source() {};
		virtual void sample() = 0;
};

class ColdGhostCell : public Source {
	public:
		ColdGhostCell(Species &species, Domain &domain, const Vector3d &x1,
				const Vector3d &x2, const Vector3d &v_drift, double n);

		virtual void sample() override;

		Species &species;
		Domain &domain;

		Vector3d x1, x2, dx, v_drift;
		double A = 1, n, n_real, V1;
};

class WarmGhostCell : public ColdGhostCell {
	public:
		WarmGhostCell(Species &species, Domain &domain, const Vector3d &x1,
				const Vector3d &x2, const Vector3d &v_drift, double n, double T);

		void sample() override;

		double T;
};

class ColdBeam : public Source {
	public:
		ColdBeam(Species &species, Domain &domain, const Vector3d &x1,
				const Vector3d &x2, const Vector3d &v_drift, double n);

		virtual void sample() override;

		Species &species;
		Domain &domain;

		Vector3d x1, x2, dx, v_drift, normal, tangent1, tangent2;
		double A = 1, n, n_real;
};

class WarmBeam : public ColdBeam {
	public:
		WarmBeam(Species &species, Domain &domain, const Vector3d &x1,
				const Vector3d &x2, const Vector3d &v_drift, double n, double T);

		void sample() override;

		double T, v_th, a;
		double V1, V2, V3; /* V1: velocity normal to the wall
							  V2,V3: velocities tangential to the wall */
};

#endif
