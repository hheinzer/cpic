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

	double nu0 = pow(QE, 4)*n*ln_Lambda
			/(8*sqrt(2)*PI*pow(EPS0, 2)*sqrt(me)*pow(K*Te, 1.5));
	double nueq = 8.0/(3*sqrt(PI))*me/mi*pow(1 + me/mi*Ti/Te, -1.5)*nu0;

	double dT0 = Te - Ti;

	double t_end = 8e-9;

	cout << "lambda_D = " << lambda_D << endl
		 << "ln_Lambda = " << ln_Lambda << endl
		 << "nu0 = " << nu0 << endl
		 << "nueq = " << nueq/nu0 << endl
		 << "Te = " << Te << endl
		 << "Ti = " << Ti << endl
		 << "Tinf = " << Tinf << endl
		 << "me = " << me << endl
		 << "mi = " << mi << endl;

	/*
	 * lambda_D = 7.88489e-07
	 * ln_Lambda = 6.82875
	 * nu0 = 1.4361e+09
	 * nueq = 4.52677e+08
	 * Te = 17406.8
	 * Ti = 8703.39
	 * Tinf = 13055.1
	 * me = 9.10938e-31
	 * mi = 3.64375e-30
	 * */

	string fname = "test/simulation/nanbu_analytic.csv";
	ofstream out(fname);
	if (!out.is_open()) {
		cerr << "Could not open '" << fname << "'" << endl;
		exit(EXIT_FAILURE);
	}

	out << "time,Te,Ti\n";

	for (double t = 0; t <= t_end; t += t_end/100.0) {
		double dT = dT0*exp(-2*nueq*t);
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
	species.push_back(Species("e-", me, -QE, 2e6, domain));
	species.push_back(Species("I+", mi,  QE, 2e6, domain));

	species[0].add_warm_box(x_min, x_max, n, {0, 0, 0}, Te);
	species[1].add_warm_box(x_min, x_max, n, {0, 0, 0}, Ti);

	vector<unique_ptr<Interaction>> interactions;
	interactions.push_back(make_unique<DSMC_Nanbu>(domain, species, Te, n));

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
		}
	}
}
