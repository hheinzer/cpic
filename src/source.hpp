#ifndef SOURCE_HPP
#define SOURCE_HPP

#include <iostream>
#include <Eigen/Dense>
#include "domain.hpp"
#include "species.hpp"

class Source {
	public:
		virtual ~Source() {};
		virtual void sample() = 0;
};

class ColdGhostCell : public Source {
	public:
		using Vector3d = Eigen::Vector3d;

		ColdGhostCell(Species &species, Domain &domain, const Vector3d &x1,
				const Vector3d &x2, const Vector3d &v_drift, double n);

		virtual void sample() override;

		Species &species;
		Domain &domain;

		Vector3d x1, x2, dx, v_drift;
		double V = 1, n, n_sim;
};

class WarmGhostCell : public ColdGhostCell {
	public:
		using Vector3d = Eigen::Vector3d;

		WarmGhostCell(Species &species, Domain &domain, const Vector3d &x1,
				const Vector3d &x2, const Vector3d &v_drift, double n, double T) :
			ColdGhostCell(species, domain, x1, x2, v_drift, n), T{T} {}

		void sample() override;

		double T;
};

#endif
