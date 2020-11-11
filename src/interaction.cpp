#include "interaction.hpp"
#include "const.hpp"

using namespace std;
using namespace Eigen;
using namespace Const;

DSMC_Bird::DSMC_Bird(Domain &domain, Species &species) :
	domain{domain}, species{species}
{
	n_cells = domain.n_cells;
	V = domain.get_del_x().prod();

	w_mp = species.w_mp0;

	m = species.m;
	mr = m*m/(m + m);

	string name = species.name;
	name.erase(remove_if(name.begin(), name.end(),
				[](char c){return !isalpha(c);}), name.end());
	d_ref = d_ref_map.at(name);
	T_ref = T_ref_map.at(name);
	omega = omega_map.at(name);
	df    = d_ref*sqrt(pow(2*K*T_ref/mr, omega - 0.5)/tgamma(2.5 - omega));
}

void DSMC_Bird::apply(double dt)
{
	vector<vector<Particle *>> pic(n_cells);
	for(Particle &p : species.particles) {
		int c = domain.x_to_c(p.x);
		pic[c].push_back(&p);
	}

	int n_collisions = 0;
	double sigma_vr_max_tmp = 0;

	for (int c = 0; c < n_cells; ++c) {
		int N_p = pic[c].size();
		if (N_p < 2) continue;

		/* Bird's No Time Counter */
		int N_g = (int)(0.5*N_p*N_p*w_mp*sigma_vr_max*dt/V + rng());

		for (int g = 0; g < N_g; ++g) {
			Particle *p1, *p2;

			p1 = pic[c][(int)(N_p*rng())];
			do {
				p2 = pic[c][(int)(N_p*rng())];
			} while(p2 == p1);

			double vr_mag = (p1->v - p2->v).norm();
			double sigma_vr = sigma(vr_mag)*vr_mag;

			if (sigma_vr > sigma_vr_max_tmp)
				sigma_vr_max_tmp = sigma_vr;

			double P = sigma_vr/sigma_vr_max;

			if (P > rng()) {
				++n_collisions;
				collide(p1->v, p2->v, m, m);
			}
		}
	}

	if (n_collisions > 0)
		sigma_vr_max = sigma_vr_max_tmp;
}

double DSMC_Bird::sigma(double vr_mag) const
{
	/* Bird's Variable Hard Sphere */
	double d = df*pow(1/(vr_mag*vr_mag), omega/2 - 0.25);
	return PI*d*d;
}

void DSMC_Bird::collide(Vector3d &v1, Vector3d &v2, double m1, double m2) const
{
	Vector3d vm = (m1*v1 + m2*v2)/(m1 + m2);

	double vr_mag = (v2 - v1).norm();

	/* isotropic scattering angle */
	double cos_xi = 2*rng() - 1;
	double sin_xi = sqrt(1- cos_xi*cos_xi);
	double eps = 2*PI*rng();

	Vector3d vr = {
		vr_mag*cos_xi,
		vr_mag*sin_xi*cos(eps),
		vr_mag*sin_xi*sin(eps)
	};

	v1 = vm + m2/(m1 + m2)*vr;
	v2 = vm - m1/(m1 + m2)*vr;
}
