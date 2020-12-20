#include "source.hpp"

using namespace std;
using namespace Eigen;
using namespace Const;

ColdGhostCell::ColdGhostCell(Species &species, Domain &domain, const Vector3d &x1,
		const Vector3d &x2, const Vector3d &v_drift, double n) :
	species{species}, domain{domain}, x1{x1}, x2{x2}, v_drift{v_drift}, n{n}
{
	dx = x2 - x1;

	/* make sure that x1 and x2 form a plane, not a volume */
	assert(dx.minCoeff() == 0);

	for(int dir : {X, Y, Z})
		if (abs(dx(dir)) > 0)
			A *= abs(dx(dir));

	assert(A > 0);

	/* calculate surface normal vector */
	int i_n; dx.minCoeff(&i_n);
	Vector3d normal;
	if (x2(i_n) == domain.get_x_min()(i_n)) {
		normal = Vector3d::Unit(i_n);
	} else {
		normal = -Vector3d::Unit(i_n);
	}

	/* expand the surface to a volume */
	this->x1 -= domain.get_del_x()(i_n)*normal;
	this->dx += domain.get_del_x()(i_n)*normal;

	V1 = v_drift.transpose()*normal;
	n_real = n*V1*A*domain.get_time_step();
}

void ColdGhostCell::sample()
{
	int n_sim = (int)(n_real/species.w_mp0 + rng());

	for (int p = 0; p < n_sim; ++p) {
		Vector3d v, x;

		do {
			v = v_drift;
			x = x1 + rng(3).cwiseProduct(dx) + v*domain.get_time_step();
		} while(!domain.is_inside(x));

		species.add_particle(x, v, domain.get_time_step());
	}
}


WarmGhostCell::WarmGhostCell(Species &species, Domain &domain, const Vector3d &x1,
		const Vector3d &x2, const Vector3d &v_drift, double n, double T) :
	ColdGhostCell(species, domain, x1, x2, v_drift, n), T{T}
{
	double v_th = sqrt(2*K*T/species.m);

	/* speed ratio */
	double a = V1/v_th;

	/* make sure it is an inflow */
	assert(a >= 0);

	/* overwrite n_real of ColdGhostCell */
	n_real = n*A*domain.get_time_step()*v_th/(2*sqrt(PI))
		*(exp(-a*a) + a*sqrt(PI)*(1 + erf(a)));
}

void WarmGhostCell::sample()
{
	int n_sim = (int)(n_real/species.w_mp0 + rng());

	for (int p = 0; p < n_sim; ++p) {
		Vector3d v, x;

		do {
			v = v_drift + species.get_maxwellian_velocity(T);
			x = x1 + rng(3).cwiseProduct(dx) + v*domain.get_time_step();
		} while(!domain.is_inside(x));

		species.add_particle(x, v, domain.get_time_step());
	}
}


ColdBeam::ColdBeam(Species &species, Domain &domain, const Vector3d &x1,
		const Vector3d &x2, const Vector3d &v_drift, double n) :
	species{species}, domain{domain}, x1{x1}, x2{x2}, v_drift{v_drift}, n{n}
{
	dx = x2 - x1;

	/* make sure that x1 and x2 form a plane, not a volume */
	assert(dx.minCoeff() == 0);

	for(int dir : {X, Y, Z})
		if (abs(dx(dir)) > 0)
			A *= abs(dx(dir));

	assert(A > 0);

	/* calculate surface normal vector */
	int i_n; dx.minCoeff(&i_n);
	if (x2(i_n) == domain.get_x_min()(i_n)) {
		normal = Vector3d::Unit(i_n);
	} else {
		normal = -Vector3d::Unit(i_n);
	}
	if (normal.cross(Vector3d::UnitX()).norm() != 0) {
		tangent1 = normal.cross(normal + Vector3d::UnitX()).normalized();
	} else {
		tangent1 = normal.cross(normal + Vector3d::UnitY()).normalized();
	}
	tangent2 = normal.cross(tangent1);

	V1 = v_drift.transpose()*normal;
	V2 = v_drift.transpose()*tangent1;
	V3 = v_drift.transpose()*tangent2;

	n_real = n*V1*A*domain.get_time_step();
}

void ColdBeam::sample()
{
	int n_sim = (int)(n_real/species.w_mp0 + rng());

	for (int p = 0; p < n_sim; ++p) {
		Vector3d v = v_drift;
		Vector3d x = x1 + rng(3).cwiseProduct(dx);
		species.add_particle(x, v, rng()*domain.get_time_step());
	}
}


WarmBeam::WarmBeam(Species &species, Domain &domain, const Vector3d &x1,
		const Vector3d &x2, const Vector3d &v_drift, double n, double T) :
	ColdBeam(species, domain, x1, x2, v_drift, n), T{T}
{
	v_th = sqrt(2*K*T/species.m);

	/* speed ratio */
	a = V1/v_th;

	/* make sure it is an inflow */
	assert(a >= 0);

	/* overwrite n_real of ColdBeam */
	n_real = n*A*domain.get_time_step()*v_th/(2*sqrt(PI))
		*(exp(-a*a) + a*sqrt(PI)*(1 + erf(a)));
}

void WarmBeam::sample()
{
	/* Garcias inflow boundary condition */
	int n_sim = (int)(n_real/species.w_mp0 + rng());

	for (int p = 0; p < n_sim; ++p) {

		double z_star;
		if (V1 > 0) {
			do {
				if (1/(2*a*sqrt(PI) + 1) > rng()) {
					z_star = -sqrt(-log(rng()));
				} else {
					z_star = 1/sqrt(2)*rng.normal();
				}
			} while((a - z_star)/a <= rng());
		} else {
			z_star = -sqrt(-log(rng()));
		}

		double v1 = V1 - z_star*v_th;
		double v2 = V2 + sqrt(0.5)*v_th*rng.normal();
		double v3 = V3 + sqrt(0.5)*v_th*rng.normal();

		Vector3d v = v1*normal + v2*tangent1 + v3*tangent2;
		Vector3d x = x1 + rng(3).cwiseProduct(dx);
		species.add_particle(x, v, rng()*domain.get_time_step());
	}
}
