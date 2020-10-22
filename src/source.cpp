#include "source.hpp"

using namespace std;
using namespace Eigen;
using namespace Const;

ColdGhostCell::ColdGhostCell(Species &species, Domain &domain, const Vector3d &x1,
		const Vector3d &x2, const Vector3d &v_drift, double n) :
	species{species}, domain{domain}, x1{x1}, x2{x2}, v_drift{v_drift}, n{n}
{
	dx = x2 - x1;

	for(int dir : {X, Y, Z})
		V *= abs(dx(dir));

	assert(V > 0);

	n_sim = n*V/species.w_mp0;
}

void ColdGhostCell::sample()
{
	for (int p = 0; p < n_sim; ++p) {
		Vector3d v = v_drift;
		Vector3d x = x1 + rng(3).cwiseProduct(dx) + v*domain.get_time_step();
		if (domain.is_inside(x))
			species.add_particle(x, v);
	}
}

void WarmGhostCell::sample()
{
	for (int p = 0; p < n_sim; ++p) {
		Vector3d v = v_drift + species.get_maxwellian_velocity(T);
		Vector3d x = x1 + rng(3).cwiseProduct(dx) + v*domain.get_time_step();
		if (domain.is_inside(x))
			species.add_particle(x, v);
	}
}
