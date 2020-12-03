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

const double n    = 1e20;			/* [1/m^3] */
const double Te   = 1.5*EvToK;		/* [K] */
const double Ti   = Te/2.0; 		/* [K] */
const double Tinf = (Te + Ti)/2;	/* [K] */
const double me   = ME;				/* [kg] */
const double mi   = 4*ME;			/* [kg] */

void save_analytical_solution()
{
	double lambda_D = sqrt(EPS0*K*Tinf/(n*QE*QE));
	double ln_Lambda = log(lambda_D*2*PI*EPS0*3*K*Tinf/(QE*QE));
	double mu = me*mi/(me + mi);
	double tau0 = 1/(n*pow(QE, 4)*ln_Lambda
			/(4*PI*pow(EPS0, 2)*sqrt(mu)*pow(K*Tinf, 1.5)));

	double t_hat_end = 50;
	double dt_hat = t_hat_end/100.0;

	double dT0 = Te - Ti;

	string fname = "test/simulation/nanbu_analytic.csv";
	ofstream out(fname);
	if (!out.is_open()) {
		cerr << "Could not open '" << fname << "'" << endl;
		exit(EXIT_FAILURE);
	}

	out << "time,Te,Ti\n";

	for (double t_hat = 0; t_hat <= t_hat_end; t_hat += dt_hat) {
		double t = t_hat*tau0;
		double dT = dT0*exp(-2/(5*sqrt(2*PI))*t_hat);

		out << t << ',' << Tinf + dT/2.0 << ',' << Tinf - dT/2.0 << endl;
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
	species.push_back(Species("e-", me, -QE, 1e6, domain));
	species.push_back(Species("I+", mi,  QE, 1e6, domain));

	species[0].add_warm_box(x_min, x_max, n, {0, 0, 0}, Te);
	species[1].add_warm_box(x_min, x_max, n, {0, 0, 0}, Ti);

	vector<unique_ptr<Interaction>> interactions;
	interactions.push_back(make_unique<DSMC_Nanbu>(domain, species, Te, n));

	//Solver solver(domain, 10000, 1e-4);

	domain.check_formulation(n, Te);

	while (domain.advance_time()) {
		//domain.calc_charge_density(species);
		//solver.calc_potential();
		//solver.calc_electric_field();

		for(Species &sp : species) {
			sp.push_particles_leapfrog();
			//sp.remove_dead_particles();
			sp.calc_number_density();
			sp.sample_moments();
			sp.calc_gas_properties();
		}

		domain.calc_total_temperature(species);

		for(auto &interaction : interactions)
			interaction->apply(domain.get_time_step());

		if (domain.get_iter()%10 == 0 || domain.is_last_iter()) {
			//for(Species &sp : species) {
			//	sp.calc_macroparticle_count();
			//}
			//domain.calc_coulomb_log(species[0]);

			domain.print_info(species);
			domain.write_statistics(species);
			//domain.save_fields(species);
			//domain.save_particles(species, 1000);
			//domain.save_velocity_histogram(species);
		}
	}
}
