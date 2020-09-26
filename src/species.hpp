#ifndef SPECIES_HPP
#define SPECIES_HPP

#include <string>
#include <iostream>
#include <Eigen/Dense>
#include "const.hpp"
#include "domain.hpp"
#include "random.hpp"

struct Particle {
	using Vector3d = Eigen::Vector3d;

	Particle(const Vector3d &x, const Vector3d &v, double dt, double w_mp) :
		x{x}, v{v}, dt{dt}, w_mp{w_mp} {}

	Vector3d x;		/* [m] particle position */
	Vector3d v;		/* [m/s] particle velocity */
	double dt;		/* [s] particle time step */
	double w_mp;	/* [-] macro particle weight */
};

class Species {
	public:
		using Vector3d = Eigen::Vector3d;
		using VectorXd = Eigen::VectorXd;
		using MatrixXd = Eigen::MatrixXd;

		Species(std::string name, double m_s, double rho_s, double w_mp0, Domain &domain);

		int get_sim_count() const {return (int)particles.size();}

		double get_real_count() const;

		Vector3d get_momentum() const;

		double get_kinetic_energy() const;

		double get_maxwellian_velocity_magnitude(double T) const;

		Vector3d get_maxwellian_velocity(double T) const;

		void add_particle(const Vector3d &x, const Vector3d &v);

		void add_particle(const Vector3d &x, const Vector3d &v, double dt, double w_mp);

		void add_cold_box(const Vector3d &x1, const Vector3d &x2, double n,
				const Vector3d &v_drift);

		void add_warm_box(const Vector3d &x1, const Vector3d &x2, double n,
				const Vector3d &v_drift, double T);

		void push_particles_leapfrog();

		void remove_dead_particles();

		void calc_number_density();

		void clear_moments();

		void sample_moments();

		void calc_gas_properties();

		void calc_macroparticle_count();

		void update_mean();

		const std::string name;
		const double m_s;	/* [kg] species mass */
		const double rho_s;	/* [C] species charge */
		const double w_mp0;	/* [-] default macro particle weight */

		std::vector<Particle> particles;
		VectorXd n;			/* [1/m^3] number density */
		VectorXd n_mean;	/* [1/m^3] time averaged number density */
		VectorXd T;
		VectorXd mp_count;
		MatrixXd v_stream;

	private:
		int n_samples = 0;

		VectorXd n_sum, nuu_sum, nvv_sum, nww_sum;
		MatrixXd nv_sum;

		Domain &domain;
};

#endif
