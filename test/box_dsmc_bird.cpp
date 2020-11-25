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

int main()
{
	Vector3d x_min, x_max, x_mid;
	x_min << -0.1, -0.1, -0.1;
	x_max <<  0.1,  0.1,  0.1;

	Domain domain("test/simulation/box_dsmc_bird", 21, 21, 21);
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
	species.push_back(Species("O", 16*AMU, 0, 1e13, domain));

	const double n = 1e20;
	const double T = 1;

	species[0].add_warm_box(x_min, x_max, n/2, { 100, 0, 0}, T);
	species[0].add_warm_box(x_min, x_max, n/2, {-100, 0, 0}, T);

	vector<unique_ptr<Interaction>> interactions;
	interactions.push_back(make_unique<DSMC_Bird>(domain, species[0]));

	domain.check_formulation(n, T);

	while (domain.advance_time()) {
		for(auto &interaction : interactions)
			interaction->apply(domain.get_time_step());

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
			domain.save_velocity_histogram(species);
		}
	}
}
