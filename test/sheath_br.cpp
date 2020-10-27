#include <vector>
#include <Eigen/Dense>
#include "const.hpp"
#include "domain.hpp"
#include "species.hpp"
#include "source.hpp"
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
	Vector3d x_min = {0.00, -0.00075, -0.00075};
	Vector3d x_max = {0.03,  0.00075,  0.00075};

	Domain domain("test/simulation/sheath_br", 21, 2, 2);
	domain.set_dimensions(x_min, x_max);
	domain.set_time_step(1e-8);
	domain.set_iter_max(2000);

	domain.set_bc_at(Xmin, BC(PBC::Open,     FBC::Dirichlet));
	domain.set_bc_at(Xmax, BC(PBC::Open,     FBC::Dirichlet, -0.18011));
	domain.set_bc_at(Ymin, BC(PBC::Periodic, FBC::Periodic));
	domain.set_bc_at(Ymax, BC(PBC::Periodic, FBC::Periodic));
	domain.set_bc_at(Zmin, BC(PBC::Periodic, FBC::Periodic));
	domain.set_bc_at(Zmax, BC(PBC::Periodic, FBC::Periodic));

	vector<Species> species;
	species.push_back(Species("O+", 16*AMU, QE, 10, domain));

	const double n = 1e12;

	vector<unique_ptr<Source>> sources;
	Vector3d x1 = {0.0, -0.00075, -0.00075};
	Vector3d x2 = {0.0,  0.00075,  0.00075};
	Vector3d vi = {11492.19, 0, 0};
	double   T  = 1000;
	sources.push_back(make_unique<WarmBeam>(species[0], domain, x1, x2, vi, n, T));

	Solver solver(domain, 1000, 1e-4);
	solver.set_reference_values(0, T*KToEv, n);

	save_analytical_solution();

	while (domain.advance_time()) {
		solver.calc_potential_BR();
		solver.calc_electric_field();

		for (auto &source : sources)
			source->sample();

		for(Species &sp : species) {
			sp.push_particles_leapfrog();
			sp.remove_dead_particles();
			sp.calc_number_density();
		}

		domain.calc_charge_density(species);

		if (domain.steady_state(species, 5)) {
			for(Species &sp : species)
				sp.update_mean();
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
			//domain.save_particles(species, 1000);
			//domain.save_velocity_histogram(species);
		}
	}
}
