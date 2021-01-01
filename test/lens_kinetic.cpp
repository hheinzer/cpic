#include <vector>
#include <Eigen/Dense>
#include "const.hpp"
#include "domain.hpp"
#include "species.hpp"
#include "source.hpp"
#include "solver.hpp"
#include "interaction.hpp"

using namespace std;
using namespace Const;
using namespace Eigen;
using PBC = ParticleBCtype;
using FBC = FieldBCtype;

int main()
{
	Vector3d x_min = {0.0, -0.05, -0.05};
	Vector3d x_max = {0.3,  0.05,  0.05};

	Domain domain("test/simulation/lens_kinetic", 61, 21, 21);
	domain.set_dimensions(x_min, x_max);
	domain.set_time_step(1e-9);
	domain.set_iter_max(30000);

	double phi_l = -100; /* [V] */

	domain.set_bc_at(Xmin, BC(PBC::Open,     FBC::Neumann));
	domain.set_bc_at(Xmax, BC(PBC::Open,     FBC::Neumann));
	domain.set_bc_at(Ymin, BC(PBC::Specular, FBC::Dirichlet));
	domain.set_bc_at(Ymax, BC(PBC::Specular, FBC::Dirichlet));
	domain.set_bc_at(Zmin, BC(PBC::Specular, FBC::Dirichlet));
	domain.set_bc_at(Zmax, BC(PBC::Specular, FBC::Dirichlet));

	auto lense = [](double x, double, double){
		return (0.1 <= x ? (x <= 0.2 ? true : false) : false); };
	domain.set_bc_at(Ymin, BC(PBC::Specular, FBC::Dirichlet, phi_l, lense));
	domain.set_bc_at(Ymax, BC(PBC::Specular, FBC::Dirichlet, phi_l, lense));
	domain.set_bc_at(Zmin, BC(PBC::Specular, FBC::Dirichlet, phi_l, lense));
	domain.set_bc_at(Zmax, BC(PBC::Specular, FBC::Dirichlet, phi_l, lense));

	vector<Species> species;
	species.push_back(Species("Xe+", 54*AMU,  QE, 1e3, domain));
	species.push_back(Species("e-",      ME, -QE, 1e3, domain));

	const double n = 1e11;

	vector<unique_ptr<Source>> sources;
	Vector3d x1 = {0.0, -0.02, -0.02};
	Vector3d x2 = {0.0,  0.02,  0.02};
	Vector3d v  = {1e4, 0, 0};
	double   T  = 1000;
	sources.push_back(make_unique<WarmBeam>(species[0], domain, x1, x2, v, n, T));
	sources.push_back(make_unique<WarmBeam>(species[1], domain, x1, x2, v, n, T));

	Solver solver(domain, 30000, 1e-4);

	domain.check_formulation(n, T);

	while (domain.advance_time()) {
		domain.calc_charge_density(species);
		solver.calc_potential();
		solver.calc_electric_field();

		for (auto &source : sources)
			source->sample();

		for (Species &sp : species) {
			sp.push_particles_leapfrog();
			sp.remove_dead_particles();
			sp.calc_number_density();
		}

		if (!domain.averaing_time() && domain.steady_state(species, 1000, 0.01)) {
			domain.start_averaging_time();
			for(Species &sp : species)
				sp.start_time_averaging(5000);
		}

		if (domain.get_iter()%1000 == 0 || domain.is_last_iter()) {
			for(Species &sp : species) {
				sp.sample_moments();
				sp.calc_gas_properties();
				sp.calc_macroparticle_count();
			}
			domain.calc_coulomb_log(T, n);
			domain.calc_total_temperature(species);

			domain.print_info(species);
			domain.write_statistics(species);
			domain.save_fields(species);
			domain.save_particles(species, 1000);
			domain.save_velocity_histogram(species);
		}
	}
}
