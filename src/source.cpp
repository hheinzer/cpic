#include "source.hpp"

using namespace std;
using namespace Eigen;
using namespace Const;

ColdBeam::ColdBeam(Species &species, Domain &domain, const Vector3d &x1,
		const Vector3d &x2, const Vector3d &v_drift, double n) :
	species{species}, domain{domain}, x1{x1}, x2{x2}, v_drift{v_drift}, n{n}
{
	dx = x2 - x1 - v_drift.normalized().cwiseProduct(domain.get_del_x());
	l_E = domain.x_to_l((x2 + x1)/2);
	V = abs(dx.prod());
	n_sim = n*V/species.w_mp0;
}

void ColdBeam::sample()
{
	Vector3d E_l = domain.gather(domain.E, l_E);
	Vector3d dv = E_l*(domain.get_time_step()*species.q/species.m);

	for (int p = 0; p < n_sim; ++p) {
		Vector3d v = v_drift + dv;
		Vector3d x = x1 + rng(3).cwiseProduct(dx) + v*domain.get_time_step();
		if (domain.is_inside(x))
			species.add_particle(x, v);
	}
}

void WarmBeam::sample()
{
	Vector3d E_l = domain.gather(domain.E, l_E);
	Vector3d dv = E_l*(domain.get_time_step()*species.q/species.m);

	for (int p = 0; p < n_sim; ++p) {
		Vector3d v = v_drift + dv + species.get_maxwellian_velocity(T);
		Vector3d x = x1 + rng(3).cwiseProduct(dx) + v*domain.get_time_step();
		if (domain.is_inside(x))
			species.add_particle(x, v);
	}
}
