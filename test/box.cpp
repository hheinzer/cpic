#include <vector>
#include <Eigen/Dense>
#include "const.hpp"
#include "domain.hpp"
#include "species.hpp"
#include "solver.hpp"

using namespace std;
using namespace Const;
using namespace Eigen;
using PBC = ParticleBCtype;
using FBC = FieldBCtype;

int main()
{
	Vector3d x_min, x_max, x_mid;
	x_min << -0.1, -0.1, -0.1;
	x_max <<  0.1,  0.1,  0.1;
	x_mid << (x_max + x_min)/2;

	Domain domain("test/box/box", 21, 21, 21);
	domain.set_dimensions(x_min, x_max);
	domain.set_time_step(2e-10);
	domain.set_iter_max(10000);

	domain.set_bc_at(Xmin, BC(PBC::Symmetric, FBC::Dirichlet));
	domain.set_bc_at(Xmax, BC(PBC::Symmetric, FBC::Dirichlet));
	domain.set_bc_at(Ymin, BC(PBC::Symmetric, FBC::Dirichlet));
	domain.set_bc_at(Ymax, BC(PBC::Symmetric, FBC::Dirichlet));
	domain.set_bc_at(Zmin, BC(PBC::Symmetric, FBC::Dirichlet));
	domain.set_bc_at(Zmax, BC(PBC::Symmetric, FBC::Dirichlet));

	vector<Species> species;
	species.push_back(Species("O+", 16*AMU,  QE, 10000, domain));
	species.push_back(Species("e-",     ME, -QE, 10000, domain));

	const double n = 1e11;

	species[0].add_cold_box(x_min, x_max, n, {0, 0, 0});
	species[1].add_cold_box(x_min, x_mid, n, {0, 0, 0});

	for(Species &sp : species)
		sp.calc_number_density();

	Solver solver(domain, 10000, 1e-4);

	while (domain.advance_time()) {
		domain.calc_charge_density(species);
		solver.calc_potential();
		solver.calc_electric_field();

		for(Species &sp : species) {
			sp.push_particles_leapfrog();
			sp.calc_number_density();
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
		}
	}
}
