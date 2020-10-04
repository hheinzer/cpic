#include <vector>
#include <Eigen/Dense>
#include "const.hpp"
#include "interaction.hpp"
#include "domain.hpp"
#include "species.hpp"
#include "source.hpp"
#include "object.hpp"
#include "solver.hpp"

using namespace std;
using namespace Const;
using namespace Eigen;
using PBC = ParticleBCtype;
using FBC = FieldBCtype;

int main()
{
	Vector3d x_min, x_max;
	x_min << 0.0, -0.1, -0.1;
	x_max << 0.4,  0.1,  0.1;

	Domain domain("test/sphere/sphere", 41, 21, 21);
	domain.set_dimensions(x_min, x_max);
	domain.set_time_step(1e-7);
	domain.set_iter_max(1000);

	double T = 1000;

	domain.set_bc_at(Xmin, BC(PBC::Open, FBC::Dirichlet));
	domain.set_bc_at(Xmax, BC(PBC::Open, FBC::Neumann));
	domain.set_bc_at(Ymin, BC(PBC::Open, FBC::Neumann));
	domain.set_bc_at(Ymax, BC(PBC::Open, FBC::Neumann));
	domain.set_bc_at(Zmin, BC(PBC::Open, FBC::Neumann));
	domain.set_bc_at(Zmax, BC(PBC::Open, FBC::Neumann));

	vector<unique_ptr<Object>> objects;
	Vector3d x_o = {0.15, 0, 0};
	objects.push_back(make_unique<Object>(domain, "test/sphere.msh", x_o, -100, T));

	domain.save_object_mesh(objects);

	exit(0);

	vector<Species> species;
	species.push_back(Species("O+", 16*AMU, QE, 1e2, domain));

	const double n = 1e10;

	vector<unique_ptr<Source>> sources;
	Vector3d x1 = {0.0, -0.1, -0.1};
	Vector3d x2 = {0.0,  0.1,  0.1};
	Vector3d  v = {7000, 0,  0};
	sources.push_back(make_unique<WarmBeam>(species[0], domain, x1, x2, v, n, T));

	//vector<unique_ptr<Interaction>> interactions;
	//interactions.push_back(make_unique<DSMC>(domain, species[0]));

	Solver solver(domain, 20000, 1e-4);
	solver.set_reference_values(0, 1.5, n);

	while (domain.advance_time()) {
		solver.calc_potential_BR();
		solver.calc_electric_field();

		for(auto &source : sources)
			source->sample();

		//for(auto &interaction : interactions)
		//	interaction->apply(domain.get_time_step());

		for(Species &sp : species) {
			sp.push_particles_leapfrog();
			sp.remove_dead_particles();
			sp.calc_number_density();
		}

		domain.calc_charge_density(species);

		if (domain.steady_state(species)) {
			for(Species &sp : species) {
				sp.update_mean();
				sp.sample_moments();
			}
		}

		if (domain.get_iter()%10 == 0 || domain.is_last_iter()) {
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
