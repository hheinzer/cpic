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
				const Vector3d &x2, const Vector3d &v0, double n);

		void sample() override;

	private:
		Species &species;
		Domain &domain;

		Vector3d x1, x2, dx, v0;
		double A, n;

		int random_dirs[2], fix_dir;
};

#endif
