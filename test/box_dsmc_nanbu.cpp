#include <vector>
#include <Eigen/Dense>
#include "const.hpp"
#include "interaction.hpp"
#include "domain.hpp"
#include "species.hpp"
#include "solver.hpp"

using namespace std;
using namespace Const;
using namespace Eigen;
using PBC = ParticleBCtype;
using FBC = FieldBCtype;

void save_analytical_solution()
{
	double Te    = 1000;		/* [K] */
	double v_I   = 11492.19;	/* [m/s] */
	double m_I   = 16*AMU;		/* [kg] */
	double n_ei  = 1e12;		/* [1/m^3] */
	double x_L   = -0.03;		/* [m] */
	double phi_0 = -0.1811;		/* [V] */

	int N = 1000;

	double y0 = m_I*v_I*v_I/(2*K*Te);

	double lambda_D = sqrt(EPS0*K*Te/(n_ei*QE*QE));
	double xi_0 = 0.0;
	double xi_L = x_L/lambda_D;
	double d_xi = xi_L/(N - 1);

	double chi_0 = -QE*phi_0/(K*Te);

	string fname = "test/sheath_br/sheath_analytic.csv";
	ofstream out(fname);
	if (!out.is_open()) {
		cerr << "Could not open '" << fname << "'" << endl;
		exit(EXIT_FAILURE);
	}

	out << "x,phi_analytic\n";

	/* Euler-Cauchy integration */
	for (int i = 1; i < N; ++i) {
		double xi  = xi_0  + d_xi;
		double chi = chi_0 + d_xi*sqrt(4*y0*(sqrt(1 + chi_0/y0) - 1) + 2*(exp(-chi_0) - 1));

		double x   = xi*lambda_D;
		double phi = -chi*K*Te/QE;

		out << -x << "," << phi << "\n";

		xi_0 = xi;
		chi_0 = chi;
	}

	out.close();
}

int main()
{
	Vector3d x_min, x_max, x_mid;
	x_min << -0.1, -0.1, -0.1;
	x_max <<  0.1,  0.1,  0.1;

	Domain domain("test/simulation/box_dsmc_nanbu", 21, 21, 21);
	domain.set_dimensions(x_min, x_max);
	domain.set_time_step(1e-7);
	domain.set_iter_max(10000);

	domain.set_bc_at(Xmin, BC(PBC::Symmetric, FBC::Dirichlet));
	domain.set_bc_at(Xmax, BC(PBC::Symmetric, FBC::Dirichlet));
	domain.set_bc_at(Ymin, BC(PBC::Symmetric, FBC::Dirichlet));
	domain.set_bc_at(Ymax, BC(PBC::Symmetric, FBC::Dirichlet));
	domain.set_bc_at(Zmin, BC(PBC::Symmetric, FBC::Dirichlet));
	domain.set_bc_at(Zmax, BC(PBC::Symmetric, FBC::Dirichlet));

	vector<Species> species;
	species.push_back(Species("e-", ME, -QE, 1e6, domain));

	const double n = 1e11;
	const double Te = 1.5*EvToK;
	const double Tx = 1.3*Te/(1.0/3.0*1.3+2.0/3.0);
	const double Ty = Te/(1.0/3.0*1.3+2.0/3.0);
	const double Tz = Te/(1.0/3.0*1.3+2.0/3.0);

	species[0].add_warm_box(x_min, x_max, n, {0, 0, 0}, {Tx, Ty, Tz});

	vector<unique_ptr<Interaction>> interactions;
	interactions.push_back(make_unique<DSMC_Nanbu>(domain, species[0]));

	//Solver solver(domain, 10000, 1e-4);

	domain.check_formulation(n, Te);

	while (domain.advance_time()) {
		//domain.calc_charge_density(species);
		//solver.calc_potential();
		//solver.calc_electric_field();

		for(auto &interaction : interactions)
			interaction->apply(domain.get_time_step());

		for(Species &sp : species) {
			sp.push_particles_leapfrog();
			//sp.remove_dead_particles();
			sp.calc_number_density();
		}

		if (domain.get_iter()%100 == 0 || domain.is_last_iter()) {
			for(Species &sp : species) {
				sp.sample_moments();
				sp.calc_gas_properties();
				sp.calc_macroparticle_count();
			}

			domain.print_info(species);
			domain.write_statistics(species);
			domain.save_fields(species);
			domain.save_particles(species, 1000);
			domain.save_velocity_histogram(species);
		}
	}
}
