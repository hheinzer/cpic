#include "interaction.hpp"
#include "const.hpp"

using namespace std;
using namespace Eigen;
using namespace Const;

DSMC::DSMC(Domain &domain, Species &sp1, Species &sp2) :
	domain{domain}, sp1{sp1}, sp2{sp2}
{
	n_cells = domain.n_cells;
	V = domain.get_del_x().prod();

	assert(sp1.w_mp0 == sp2.w_mp0);
	w_mp = sp1.w_mp0;

	m1 = sp1.m_s;
	m2 = sp2.m_s;
	mr = m1*m2/(m1 + m2);

	string name1 = sp1.name;
	name1.erase(remove_if(name1.begin(), name1.end(),
				[](char c){return !isalpha(c);}), name1.end());
	d_ref1 = d_ref.at(name1);
	T_ref1 = T_ref.at(name1);
	omega1 = omega.at(name1);
	df1    = d_ref1*sqrt(pow(2*K*T_ref1/mr, omega1 - 0.5)/tgamma(2.5 - omega1));

	string name2 = sp2.name;
	name2.erase(remove_if(name2.begin(), name2.end(),
				[](char c){return !isalpha(c);}), name2.end());
	d_ref2 = d_ref.at(name2);
	T_ref2 = T_ref.at(name2);
	omega2 = omega.at(name2);
	df2    = d_ref2*sqrt(pow(2*K*T_ref2/mr, omega2 - 0.5)/tgamma(2.5 - omega2));
}

void DSMC::apply(double dt)
{
	vector<vector<Particle *>> pic1(n_cells);
	for(Particle &p : sp1.particles) {
		int c = domain.x_to_c(p.x);
		pic1[c].push_back(&p);
	}

	vector<vector<Particle *>> pic2(n_cells);
	for(Particle &p : sp2.particles) {
		int c = domain.x_to_c(p.x);
		pic2[c].push_back(&p);
	}

	int n_collisions = 0;
	double sigma_vr_max_tmp = 0;

	for (int c = 0; c < n_cells; ++c) {
		int N_1 = pic1[c].size();
		int N_2 = pic2[c].size();

		if (N_1 < 1 || N_2 < 1) continue;

		/* Bird's No Time Counter */
		int N_g = (int)(0.5*N_1*N_2*w_mp*sigma_vr_max*dt/V + rng());

		for (int g = 0; g < N_g; ++g) {
			Particle *p1, *p2;
			p1 = pic1[c][(int)(N_1*rng())];

			do {
				p2 = pic2[c][(int)(N_2*rng())];
			} while(p2 == p1);

			double vr = (p1->v - p2->v).norm();
			double sigma_vr = sigma(vr)*vr;

			if (sigma_vr > sigma_vr_max_tmp)
				sigma_vr_max_tmp = sigma_vr;

			double P = sigma_vr/sigma_vr_max;

			if (P > rng()) {
				++n_collisions;
				collide(p1->v, p2->v, m1, m2);
			}
		}
	}

	if (n_collisions > 0)
		sigma_vr_max = sigma_vr_max_tmp;
}

double DSMC::sigma(double v_r) const
{
	/* Bird's Variable Hard Sphere */
	double d1 = df1*pow(1/(v_r*v_r), omega1/2 - 0.25);
	double d2 = df2*pow(1/(v_r*v_r), omega2/2 - 0.25);
	double d12 = (d1 + d2)/2;
	return PI*d12*d12;
}

void DSMC::collide(Vector3d &v1, Vector3d &v2, double m1, double m2) const
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
