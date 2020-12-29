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
	Vector3d x_min = {0.0, -0.025, -0.025};
	Vector3d x_max = {0.3,  0.025,  0.025};

	Domain domain("test/simulation/lens", 31, 11, 11);
	domain.set_dimensions(x_min, x_max);
	domain.set_time_step(2e-9);
	domain.set_iter_max(10000);

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
	sources.push_back(make_unique<WarmGhostCell>(species[0], domain, x1, x2, v, n, T));
	sources.push_back(make_unique<WarmGhostCell>(species[1], domain, x1, x2, v, n, T));

	Solver solver(domain, 10000, 1e-4);
	solver.set_reference_values(0.0, T, n);

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

		if (!domain.averaing_time() && domain.steady_state(species, 50, 0.01)) {
			domain.start_averaging_time();
			for(Species &sp : species)
				sp.start_time_averaging(100);
		}

		if (domain.get_iter()%100 == 0 || domain.is_last_iter()) {
			for(Species &sp : species) {
				sp.sample_moments();
				sp.calc_gas_properties();
				sp.calc_macroparticle_count();
			}
			domain.calc_coulomb_log(species[1]);

			domain.print_info(species);
			domain.write_statistics(species);
			domain.save_fields(species);
			domain.save_particles(species, 1000);
			domain.save_velocity_histogram(species);
		}
	}
}
