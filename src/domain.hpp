#ifndef DOMAIN_HPP
#define DOMAIN_HPP

#include <map>
#include <memory>
#include <chrono>
#include <string>
#include <vector>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <Eigen/Eigen>

enum ChartesianDirection {X, Y, Z, W};

enum BoundarySide {Xmin, Xmax, Ymin, Ymax, Zmin, Zmax};

enum class FieldBCtype {Dirichlet, Neumann, Periodic};

enum class ParticleBCtype {Specular, Open, Diffuse, Symmetric, Periodic};

class Particle;
class Species;
class Object;

class BC {
	public:
		BC() {}

		BC(ParticleBCtype pbct, FieldBCtype fbct) :
			particle_bc_type{pbct}, field_bc_type{fbct} {}

		BC(ParticleBCtype pbct, FieldBCtype fbct, double value) :
			particle_bc_type{pbct}, field_bc_type{fbct}, value{value} {}

		BC(ParticleBCtype pbct, double T, FieldBCtype fbct) :
			particle_bc_type{pbct}, field_bc_type{fbct}, T{T} {}

		BC(ParticleBCtype pbct, double T, double a_th, FieldBCtype fbct) :
			particle_bc_type{pbct}, field_bc_type{fbct}, T{T}, a_th{a_th} {}

		BC(ParticleBCtype pbct, double T, FieldBCtype fbct, double value) :
			particle_bc_type{pbct}, field_bc_type{fbct}, T{T}, value{value} {}

		BC(ParticleBCtype pbct, double T, double a_th, FieldBCtype fbct, double value) :
			particle_bc_type{pbct}, field_bc_type{fbct}, T{T}, a_th{a_th}, value{value} {}

		double get_value() const {return value;}

		void set_delta(double delta) {this->value *= delta;}

		const ParticleBCtype particle_bc_type = ParticleBCtype::Specular;
		const FieldBCtype field_bc_type = FieldBCtype::Dirichlet;

		const double T = 1000;
		const double a_th = 1;	/* thermal accomodation coefficient */

	private:
		double value = 0;
};

class Domain {
	public:
		using Vector3i = Eigen::Vector3i;
		using Vector3d = Eigen::Vector3d;
		using VectorXd = Eigen::VectorXd;
		using MatrixXd = Eigen::MatrixXd;
		using T = Eigen::Triplet<double>;

		Domain(std::string prefix, int ni, int nj, int nk);

		~Domain();

		void set_dimensions(const Vector3d &x_min, const Vector3d &x_max);

		void set_time_step(double dt) {this->dt = dt;}

		void set_iter_max(int iter_max) {this->iter_max = iter_max;}

		void set_bc_at(BoundarySide side, BC bc);

		Vector3d get_x_min() const {return x_min;}

		Vector3d get_x_max() const {return x_max;}

		Vector3d get_del_x() const {return del_x;}

		double get_time_step() const {return dt;}

		int get_iter() const {return iter;}

		double get_wtime() const;

		double get_potential_energy() const;

		bool is_last_iter() const {return iter == iter_max;}

		bool is_inside(const Vector3d &x) const;

		bool advance_time() {time += dt; ++iter; return iter <= iter_max;}

		Vector3d x_to_l(const Vector3d &x) const;

		int x_to_c(const Vector3d &x) const;

		int at(int i, int j, int k) const {return i + j*ni + k*ni*nj;}

		void scatter(VectorXd &f, const Vector3d &l, double value);

		void scatter(MatrixXd &f, const Vector3d &l, const Vector3d &value);

		Vector3d gather(const MatrixXd &f, const Vector3d &l) const;

		void calc_charge_density(std::vector<Species> &species);

		void apply_boundary_conditions(const Species &sp, const Vector3d &x_old,
				Particle &p) const;

		void eval_field_BC(BoundarySide side, VectorXd &b0, std::vector<T> &coeffs,
				int u, int v) const;

		bool is_periodic(BoundarySide side) const;

		bool steady_state(std::vector<Species> &species, int check_every, double tol);

		bool steady_state() const {return is_steady_state;}

		bool averaing_time() const {return is_averaing_time;}

		void start_averaging_time() {this->is_averaing_time = true;}

		void check_formulation(double n_e, double T_e, const std::vector<double> &z_i,
				const std::vector<double> &n_i, const std::vector<double> &T_i) const;

		void print_info(std::vector<Species> &species) const;

		void write_statistics(std::vector<Species> &species);

		void save_fields(std::vector<Species> &species) const;

		void save_particles(std::vector<Species> &species, int n_particles) const;

		void save_velocity_histogram(std::vector<Species> &species) const;

		const std::string prefix;
		const int ni, nj, nk;
		const Vector3i nn;
		const int n_nodes, n_cells;

		VectorXd V_node;	/* [m^3] node volume */
		VectorXd rho;		/* [C] charge density */
		VectorXd phi;		/* [V] electric potential */
		MatrixXd E;			/* [V/m] electric field */
		VectorXd n_e_BR;	/* [1/m^3] electron density (Boltzmann relation) */

	private:
		void calc_node_volume();

		void eval_particle_BC(const Species &sp, int side, const Vector3d &X,
				const Vector3d &x_old, Particle &p, int dim, const Vector3d &n) const;

		Vector3d get_diffuse_vector(const Vector3d &n) const;

		Vector3d x_min, x_max, del_x;

		std::map<int, std::unique_ptr<BC>> bc;

		double time = 0, dt;
		int iter = -1, iter_max;

		std::chrono::time_point<std::chrono::high_resolution_clock> wtime_start;

		bool is_steady_state = false, is_averaing_time = false;
		double prev_n_tot = 0, prev_I_tot = 0, prev_E_tot = 0;

		std::ofstream stats;

		Vector3d L;		/* [m] domain lengths */
};

#endif
