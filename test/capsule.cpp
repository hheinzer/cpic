#include <vector>
#include <Eigen/Dense>
#include "const.hpp"
#include "interaction.hpp"
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
	Vector3d x_min, x_max, x_mid;
	x_min << 0.0, -0.15, -0.15;
	x_max << 0.5,  0.15,  0.15;

	Domain domain("test/capsule/capsule", 51, 31, 31);
	domain.set_dimensions(x_min, x_max);
	domain.set_time_step(1e-7);
	domain.set_iter_max(5000);

	double T = 1000;

	domain.set_bc_at(Xmin, BC(PBC::Open,       FBC::Dirichlet));
	domain.set_bc_at(Xmax, BC(PBC::Diffuse, T, FBC::Dirichlet));
	domain.set_bc_at(Ymin, BC(PBC::Open,       FBC::Neumann));
	domain.set_bc_at(Ymax, BC(PBC::Open,       FBC::Neumann));
	domain.set_bc_at(Zmin, BC(PBC::Open,       FBC::Neumann));
	domain.set_bc_at(Zmax, BC(PBC::Open,       FBC::Neumann));

	vector<Species> species;
	species.push_back(Species("O", 16*AMU, 0, 1e13, domain));

	const double n = 1e20;

	vector<unique_ptr<Source>> sources;
	Vector3d x1 = {0.0, -0.1, -0.1};
	Vector3d x2 = {0.0,  0.1,  0.1};
	Vector3d  v = {10000, 0,  0};
	sources.push_back(make_unique<WarmBeam>(species[0], domain, x1, x2, v, n, T));

	vector<unique_ptr<Interaction>> interactions;
	interactions.push_back(make_unique<DSMC>(domain, species[0]));

	while (domain.advance_time()) {
		for(auto &source : sources)
			source->sample();

		for(auto &interaction : interactions)
			interaction->apply(domain.get_time_step());

		for(Species &sp : species) {
			sp.push_particles_leapfrog();
			sp.remove_dead_particles();
			sp.calc_number_density();
		}

		if (domain.get_iter()%50 == 0 && domain.steady_state(species)) {
			for(Species &sp : species) {
				sp.update_mean();
				sp.sample_moments();
			}
		}

		if (domain.get_iter()%100 == 0 || domain.is_last_iter()) {
			for(Species &sp : species) {
				if (!domain.steady_state()) {
					sp.clear_moments();
					sp.sample_moments();
				}
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
