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

class ColdBeam : public Source {
	public:
		using Vector3d = Eigen::Vector3d;

		ColdBeam(Species &species, Domain &domain, const Vector3d &x1,
				const Vector3d &x2, const Vector3d &v_drift, double n);

		virtual void sample() override;

		Species &species;
		Domain &domain;

		Vector3d x1, x2, dx, l_E, v_drift;
		double V, n, n_sim;
};

class WarmBeam : public ColdBeam {
	public:
		using Vector3d = Eigen::Vector3d;

		WarmBeam(Species &species, Domain &domain, const Vector3d &x1,
				const Vector3d &x2, const Vector3d &v_drift, double n, double T) :
			ColdBeam(species, domain, x1, x2, v_drift, n), T{T} {}

		void sample() override;

		double T;
};

#endif
