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


DSMC_Nanbu::DSMC_Nanbu(Domain &domain, vector<Species> &species, double T_e,
		double n_e) :
	domain{domain}, species{species}
{
	n_species = species.size();
	n_cells = domain.n_cells;

	lambda_D = sqrt(EPS0*K*T_e/(n_e*QE*QE));
}

void DSMC_Nanbu::apply(double dt)
{
	vector<vector<vector<Particle *>>> sic(n_species);
	for(int s = 0; s < n_species; ++s) {
		vector<vector<Particle *>> pic(n_cells);

		for(Particle &p : species[s].particles) {
			int c = domain.x_to_c(p.x);
			pic[c].push_back(&p);
		}

		sic[s] = pic;
	}

	/* perform like collisions */
	for (int s = 0; s < n_species; ++s) {
		for (int c = 0; c < n_cells; ++c) {

			/* reference to the particles of species s in the cell c */
			vector<Particle *> &pic = sic[s][c];

			/* number particles */
			int N = pic.size();

			if (N > 1) {

				/* shuffle the particles of species s in cell c */
				std::shuffle(pic.begin(), pic.end(), rng.get_gen());

				/* get total temperature in cell c */
				double T_tot = domain.T_tot(c);

				/* collide N/2 parts */
				for (int i = 0; i + 1 < N; i += 2)
					collide(pic[i]->v, pic[i + 1]->v, species[s].m, species[s].m,
							T_tot, species[s].q, species[s].q, species[s].n_mean[c], dt);


				/* handle odd particle numbers */
				if (N%2 != 0)
					collide(pic[N - 1]->v, pic[0]->v, species[s].m, species[s].m,
							T_tot, species[s].q, species[s].q, species[s].n_mean[c], dt);
			}

		}
	}

	/* perform unlike collisions */
}

void DSMC_Nanbu::collide(Vector3d &v1, Vector3d &v2, double m1, double m2,
		double T_tot, double q1, double q2, double n2, double dt) const
{
	/* relative velocity */
	Vector3d g = v1 - v2;
	double g_mag = g.norm();
	double g_perp = sqrt(g(Y)*g(Y) + g(Z)*g(Z));

	/* calculate coulomb logarithm */
	double ln_Lambda = log(lambda_D*2*PI*EPS0*3*K*T_tot/fabs(q1*q2));

	/* calculate mass ratio */
	double mu = m1*m2/(m1 + m2);

	/* calculate collision parameter */
	double s = ln_Lambda/(4*PI)*pow(q1*q2/(EPS0*mu), 2)*n2*pow(g_mag, -3)*dt;

	/* calculate cosine and sine of scattering angle */
	double cos_xi;
	if (s < 0.1) {
		cos_xi = 1 + s*log(rng());
	} else if (0.1 <= s && s < 3.0) {
		double A_inv = 0.0056958 + 0.9560202*s - 0.508139*s*s
			+ 0.47913906*s*s*s - 0.12788975*s*s*s*s + 0.02389567*s*s*s*s*s;

		double A = 1.0/A_inv;

		cos_xi = A_inv*log(exp(-A) + 2.0*rng()*sinh(A));
	} else if (3.0 <= s && s < 6.0) {
		double A = 3.0*exp(-s);

		cos_xi = 1.0/A*log(exp(-A) + 2.0*rng()*sinh(A));
	} else {
		cos_xi = 2.0*rng() - 1.0;
	}

	double sin_xi = sqrt(1 - cos_xi*cos_xi);

	/* calculate binary collision parameter */
	double eps = 2*PI*rng();
	Vector3d h = {
		g_perp*cos(eps),
		-(g(Y)*g(X)*cos(eps) + g_mag*g(Z)*sin(eps))/g_perp,
		-(g(Z)*g(X)*cos(eps) - g_mag*g(Y)*sin(eps))/g_perp
	};

	/* perfom binary collision */
	v1 -= m2/(m1 + m2)*(g*(1 - cos_xi) + h*sin_xi);
	v2 += m1/(m1 + m2)*(g*(1 - cos_xi) + h*sin_xi);
}
