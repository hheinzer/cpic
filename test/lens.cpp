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

int main()
{
	Vector3d x_min = {0.0, -0.05, -0.05};
	Vector3d x_max = {0.4,  0.05,  0.05};

	Domain domain("test/simulation/lens", 81, 21, 21);
	domain.set_dimensions(x_min, x_max);
	domain.set_time_step(1e-7);
	domain.set_iter_max(1000);

	double phi0 = 1315; /* [V] */
	double phi1 = 1300; /* [V] */
	double phi2 = -300; /* [V] */

	domain.set_bc_at(Xmin, BC(PBC::Open,     FBC::Dirichlet, phi0));
	domain.set_bc_at(Xmax, BC(PBC::Open,     FBC::Neumann));
	domain.set_bc_at(Ymin, BC(PBC::Specular, FBC::Neumann));
	domain.set_bc_at(Ymax, BC(PBC::Specular, FBC::Neumann));
	domain.set_bc_at(Zmin, BC(PBC::Specular, FBC::Neumann));
	domain.set_bc_at(Zmax, BC(PBC::Specular, FBC::Neumann));

	/* first lens */
	domain.set_bc_at(Ymin, BC(PBC::Specular, FBC::Dirichlet, phi1, [](double x, double, double){ return (0.125 <= x ? (x <= 0.175 ? true : false) : false); }));
	domain.set_bc_at(Ymax, BC(PBC::Specular, FBC::Dirichlet, phi1, [](double x, double, double){ return (0.125 <= x ? (x <= 0.175 ? true : false) : false); }));
	domain.set_bc_at(Zmin, BC(PBC::Specular, FBC::Dirichlet, phi1, [](double x, double, double){ return (0.125 <= x ? (x <= 0.175 ? true : false) : false); }));
	domain.set_bc_at(Zmax, BC(PBC::Specular, FBC::Dirichlet, phi1, [](double x, double, double){ return (0.125 <= x ? (x <= 0.175 ? true : false) : false); }));

	/* second lens */
	domain.set_bc_at(Ymin, BC(PBC::Specular, FBC::Dirichlet, phi2, [](double x, double, double){ return (0.225 <= x ? (x <= 0.275 ? true : false) : false); }));
	domain.set_bc_at(Ymax, BC(PBC::Specular, FBC::Dirichlet, phi2, [](double x, double, double){ return (0.225 <= x ? (x <= 0.275 ? true : false) : false); }));
	domain.set_bc_at(Zmin, BC(PBC::Specular, FBC::Dirichlet, phi2, [](double x, double, double){ return (0.225 <= x ? (x <= 0.275 ? true : false) : false); }));
	domain.set_bc_at(Zmax, BC(PBC::Specular, FBC::Dirichlet, phi2, [](double x, double, double){ return (0.225 <= x ? (x <= 0.275 ? true : false) : false); }));

	vector<Species> species;
	species.push_back(Species("Xe+", 54*AMU,  QE, 1000, domain));
	//species.push_back(Species("e-",     ME, -QE, 100000, domain));

	const double n = 1e11;

	vector<unique_ptr<Source>> sources;
	Vector3d x1 = {0.0, -0.05, -0.05};
	Vector3d x2 = {0.0,  0.05,  0.05};
	Vector3d v  = {2500, 0, 0};
	double   Ti = 450;
	double   Te = 3.5*EvToK;
	sources.push_back(make_unique<WarmBeam>(species[0], domain, x1, x2, v, n, Ti));
	//sources.push_back(make_unique<WarmGhostCell>(species[1], domain, x1, x2, v, n, Te));

	Solver solver(domain, 10000, 1e-4);
	solver.set_reference_values(phi0, Te, n);

	//domain.check_formulation(n, T, {1}, {n}, {T});

	while (domain.advance_time()) {
		solver.calc_potential_BR();
		solver.calc_electric_field();

		for (auto &source : sources)
			source->sample();

		for (Species &sp : species) {
			sp.push_particles_leapfrog();
			sp.remove_dead_particles();
			sp.calc_number_density();
		}

		domain.calc_charge_density(species);

		if (!domain.averaing_time() && domain.steady_state(species, 50, 0.01)) {
			domain.start_averaging_time();
			for(Species &sp : species)
				sp.start_time_averaging(100);
		}

		if (domain.get_iter()%50 == 0 || domain.is_last_iter()) {
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
