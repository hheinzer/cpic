#include "species.hpp"
#include "random.hpp"

using namespace std;
using namespace Eigen;
using namespace Const;

Species::Species(string name, double m_s, double rho_s, double w_mp0, Domain &domain) :
	name{name}, m_s{m_s}, rho_s{rho_s}, w_mp0{w_mp0}, domain{domain}
{
	int n_nodes = domain.n_nodes;
	int n_cells = domain.n_cells;
	n        = VectorXd::Zero(n_nodes);
	n_mean   = VectorXd::Zero(n_nodes);
	T        = VectorXd::Zero(n_nodes);
	v_stream = MatrixXd::Zero(n_nodes, 3);
	mp_count = VectorXd::Zero(n_cells);
	n_sum    = VectorXd::Zero(n_nodes);
	nv_sum   = MatrixXd::Zero(n_nodes, 3);
	nuu_sum  = VectorXd::Zero(n_nodes);
	nvv_sum  = VectorXd::Zero(n_nodes);
	nww_sum  = VectorXd::Zero(n_nodes);
}

double Species::get_real_count() const
{
	double w_mp_sum = 0;
	for(const Particle &p : particles)
		w_mp_sum += p.w_mp;
	return w_mp_sum;
}

Vector3d Species::get_momentum() const
{
	Vector3d I;
	for(const Particle &p : particles)
		I += p.w_mp*p.v;
	return m_s*I;
}

double Species::get_kinetic_energy() const
{
	double E_kin = 0;
	for(const Particle &p : particles)
		E_kin += p.w_mp*p.v.squaredNorm();
	return 0.5*m_s*E_kin;
}

void Species::add_particle(const Vector3d &x, const Vector3d &v)
{
	add_particle(x, v, w_mp0);
}

void Species::add_particle(const Vector3d &x, const Vector3d &v, double w_mp)
{
	Vector3d l = domain.x_to_l(x);
	Vector3d E_p = domain.gather(domain.E, l);
	Vector3d dv = rho_s/m_s*E_p*0.5*domain.get_time_step();
	particles.emplace_back(x, v - dv, w_mp);
}

void Species::add_particle_box(const Vector3d &x1, const Vector3d &x2, double n)
{
	double V_box = (x2 - x1).prod();
	double n_real = n*V_box;
	int n_sim = (int)(n_real/w_mp0);

	for (int p = 0; p < n_sim; ++p) {
		Vector3d x;
		do {
			x = x1.array() + rng(3).array()*(x2 - x1).array();
		} while(!domain.is_inside(x));
		add_particle(x, Vector3d::Zero(), w_mp0);
	}
}

void Species::push_particles_leapfrog()
{
	double dt = domain.get_time_step();

	for(Particle &p : particles) {
		Vector3d l = domain.x_to_l(p.x);
		Vector3d E_p = domain.gather(domain.E, l);

		p.v += E_p*(dt*rho_s/m_s);
		p.x += p.v*dt;

		while (!domain.is_inside(p.x) && p.w_mp > 0)
			domain.apply_boundary_conditions(p);
	}
}

void Species::remove_dead_particles()
{
	int n_sim = get_sim_count();
	for (int p = 0; p < n_sim; ++p) {
		if (particles[p].w_mp > 0) continue;
		particles[p] = particles[n_sim - 1];
		--n_sim;
		--p;
	}
	particles.erase(particles.begin() + n_sim, particles.end());
}

void Species::calc_number_density()
{
	n.setZero();
	for(const Particle &p : particles) {
		Vector3d l = domain.x_to_l(p.x);
		domain.scatter(n, l, p.w_mp);
	}
	n = n.array()/domain.V_node.array();
}

void Species::sample_moments()
{
	for(const Particle &p : particles) {
		Vector3d l = domain.x_to_l(p.x);
		domain.scatter(n_sum,   l, p.w_mp);
		domain.scatter(nv_sum,  l, p.w_mp*p.v);
		domain.scatter(nuu_sum, l, p.w_mp*p.v(X)*p.v(X));
		domain.scatter(nvv_sum, l, p.w_mp*p.v(Y)*p.v(Y));
		domain.scatter(nww_sum, l, p.w_mp*p.v(Z)*p.v(Z));
	}
}

void Species::calc_gas_properties()
{
	for (int u = 0; u < domain.n_nodes; ++u) {
		double count = n_sum(u);
		if (count <= 0) {
			v_stream.row(u).setZero();
			T(u) = 0;
			continue;
		}

		v_stream.row(u) = nv_sum.row(u)/count;

		double u_mean = v_stream(u, X);
		double v_mean = v_stream(u, Y);
		double w_mean = v_stream(u, Z);

		double u2_mean = nuu_sum(u)/count;
		double v2_mean = nvv_sum(u)/count;
		double w2_mean = nww_sum(u)/count;

		double uu = u2_mean - u_mean*u_mean;
		double vv = v2_mean - v_mean*v_mean;
		double ww = w2_mean - w_mean*w_mean;

		T(u) = m_s/(2*K)*(uu + vv + ww);
	}
}

void Species::calc_macroparticle_density()
{
	mp_count.setZero();
	for(const Particle &p : particles) {
		int c = domain.x_to_c(p.x);
		mp_count(c) += 1;
	}
}

void Species::clear_samples()
{
	n_sum.setZero();
	nv_sum.setZero();
	nuu_sum.setZero();
	nvv_sum.setZero();
	nww_sum.setZero();
}

void Species::update_mean()
{
	n_mean = (n.array() + n_samples*n_mean.array())/(n_samples + 1);
	++n_samples;
}

double Species::sample_thermal_vel(double T) const
{
	double v_th = sqrt(2*K*T/m_s);
	Vector3d v_M = v_th*(rng(3, M).rowwise().sum().array() - M/2.0)*sqrt(6.0/M);
	return v_M.norm();
}

Vector3d Species::sample_isothermal_vel(double T) const
{
	double theta = 2*PI*rng();
	double r = -1 + 2*rng();
	double a = sqrt(1 - r*r);
	Vector3d d(r, cos(theta)*a, sin(theta)*a);
	double v_th = sample_thermal_vel(T);
	return v_th*d;
}

Vector3d Species::sample_reflected_vel(const Vector3d &d, double v_mag1, double T) const
{
	double v_th = sample_thermal_vel(T);
	static const double a_th = 1;
	double v_mag2 = v_mag1 + a_th*(v_th - v_mag1);
	return v_mag2*d;
}
