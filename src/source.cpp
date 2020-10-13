#include "source.hpp"

using namespace std;
using namespace Eigen;
using namespace Const;

ColdBeam::ColdBeam(Species &species, Domain &domain, const Vector3d &x1,
		const Vector3d &x2, const Vector3d &v_drift, double n) :
	species{species}, domain{domain}, x1{x1}, x2{x2}, v_drift{v_drift}, n{n}
{
	dx = x2 - x1;

	for(int dir : {X, Y, Z})
		if (dx(dir) != 0)
			A *= dx(dir);
}

void ColdBeam::sample()
{
	double n_real = n*v_drift.norm()*A*domain.get_time_step();
	int n_sim = (int)(n_real/species.w_mp0 + rng());

	for (int p = 0; p < n_sim; ++p) {
		Vector3d v = v_drift;
		Vector3d x = x1 + rng(3).cwiseProduct(dx);
		species.add_particle(x, v);
	}
}

void WarmBeam::sample()
{
	double v_th = sqrt(K*T/(2*PI*species.m));
	double n_real = n*(v_drift.norm() + v_th)*A*domain.get_time_step();
	int n_sim = (int)(n_real/species.w_mp0 + rng());

	for (int p = 0; p < n_sim; ++p) {
		Vector3d v = v_drift + species.get_maxwellian_velocity(T);
		Vector3d x = x1 + rng(3).cwiseProduct(dx);
		species.add_particle(x, v);
	}
}
