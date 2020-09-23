#include <vector>
#include <Eigen/Dense>
#include "const.hpp"
#include "domain.hpp"
#include "species.hpp"
#include "solver.hpp"

using namespace std;
using namespace Const;
using namespace Eigen;

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

	domain.set_boundary_condition(Xmin, BC(Dirichlet, Open));
	domain.set_boundary_condition(Xmax, BC(Dirichlet, Open));
	domain.set_boundary_condition(Ymin, BC(Dirichlet, Open));
	domain.set_boundary_condition(Ymax, BC(Dirichlet, Open));
	domain.set_boundary_condition(Zmin, BC(Dirichlet, Open));
	domain.set_boundary_condition(Zmax, BC(Dirichlet, Open));

	vector<Species> species;
	species.push_back(Species("O+", 16*AMU,  QE, 8000, domain));
	species.push_back(Species("e-",     ME, -QE, 1000, domain));

	species[0].add_particle_box(0.5*x_min, 0.5*x_max, 1e11);
	species[1].add_particle_box(0.5*x_min,     x_mid, 1e11);

	for(Species &sp : species)
		sp.calc_number_density();
	domain.calc_charge_density(species);

	Solver solver(domain, 10000, 1e-4);
	solver.calc_potential();
	solver.calc_electric_field();

	while (domain.advance_time()) {
		for(Species &sp : species) {
			sp.push_particles_leapfrog();
			sp.remove_dead_particles();
			sp.calc_number_density();
			sp.sample_moments();
		}

		domain.calc_charge_density(species);

		solver.calc_potential();
		solver.calc_electric_field();

		if (domain.get_iter()%100 == 0 || domain.is_last_iter()) {
			domain.print_info(species);
			domain.write_statistics(species);
			domain.save_fields(species);
		}
	}
}
