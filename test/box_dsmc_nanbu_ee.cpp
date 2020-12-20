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

const double ne = 1e20;							/* [1/m^3] */
const double Te = 1.5*EvToK;					/* [K] */
const double Ty = Te/(1.0/3.0*1.3+2.0/3.0);		/* [K] */
const double Tx = 1.3*Ty;						/* [K] */

void save_analytical_solution()
{
	double lambda_D = sqrt(EPS0*K*Te/(ne*QE*QE));
	double ln_Lambda = log(lambda_D*2*PI*EPS0*3*K*Te/(QE*QE));
	double tau0 = 1/((ne*pow(QE, 4)*ln_Lambda
			/(8*PI*sqrt(2)*pow(EPS0, 2)*sqrt(ME)*pow(K*Te, 1.5))));

	double t_hat_end = 12;
	double dt_hat = 12.0/100.0;

	double dT0 = Tx - Ty;

	string fname = "test/simulation/nanbu_analytic.csv";
	ofstream out(fname);
	if (!out.is_open()) {
		cerr << "Could not open '" << fname << "'" << endl;
		exit(EXIT_FAILURE);
	}

	out << "t,Tx,Ty\n";

	for (double t_hat = 0; t_hat <= t_hat_end; t_hat += dt_hat) {
		double t = t_hat*tau0;
		double dT = dT0*exp(-8/(5*sqrt(2*PI))*t_hat);

		out << t << ',' << Te + 2.0/3.0*dT << ',' << Te - 1.0/3.0*dT  << endl;
	}

	out.close();
}

int main()
{
	save_analytical_solution();

	Vector3d x_min, x_max, x_mid;
	x_min << -0.0005, -0.0005, -0.0005;
	x_max <<  0.0005,  0.0005,  0.0005;

	Domain domain("test/simulation/box_dsmc_nanbu", 2, 2, 2);
	domain.set_dimensions(x_min, x_max);
	domain.set_time_step(1e-11);
	domain.set_iter_max(800);

	domain.set_bc_at(Xmin, BC(PBC::Periodic, FBC::Periodic));
	domain.set_bc_at(Xmax, BC(PBC::Periodic, FBC::Periodic));
	domain.set_bc_at(Ymin, BC(PBC::Periodic, FBC::Periodic));
	domain.set_bc_at(Ymax, BC(PBC::Periodic, FBC::Periodic));
	domain.set_bc_at(Zmin, BC(PBC::Periodic, FBC::Periodic));
	domain.set_bc_at(Zmax, BC(PBC::Periodic, FBC::Periodic));

	vector<Species> species;
	species.push_back(Species("e-", ME, -QE, 0.5e6, domain));

	species[0].add_warm_box(x_min, x_max, ne, {1e6, 0, 0}, {Tx, Ty, Ty});

	vector<unique_ptr<Interaction>> interactions;
	interactions.push_back(make_unique<DSMC_Nanbu>(domain, species, Te, ne));

	while (domain.advance_time()) {
		for(Species &sp : species) {
			sp.push_particles_leapfrog();
			sp.calc_number_density();
			sp.sample_moments();
			sp.calc_gas_properties();
		}

		domain.calc_total_temperature(species);

		for(auto &interaction : interactions)
			interaction->apply(domain.get_time_step());

		if (domain.get_iter()%10 == 0 || domain.is_last_iter()) {
			domain.print_info(species);
			domain.write_statistics(species);
			domain.save_velocity_histogram(species);
		}
	}
}
