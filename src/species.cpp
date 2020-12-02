#include "species.hpp"
#include "random.hpp"

using namespace std;
using namespace Eigen;
using namespace Const;

Species::Species(string name, double m, double q, double w_mp0, Domain &domain) :
	name{name}, m{m}, q{q}, w_mp0{w_mp0}, domain{domain}
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
	Vector3d I = Vector3d::Zero();
	for(const Particle &p : particles)
		I += p.w_mp*p.v;
	return m*I;
}

double Species::get_kinetic_energy() const
{
	double E_kin = 0;
	for(const Particle &p : particles)
		E_kin += p.w_mp*p.v.squaredNorm();
	return 0.5*m*E_kin;
}

Vector3d Species::get_translation_temperature() const
{
	Vector3d c2_mean = Vector3d::Zero();
	for(const Particle &p : particles)
		c2_mean.array() += p.v.array().pow(2);
	return m/(K*get_sim_count())*c2_mean;
}

double Species::get_maxwellian_velocity_magnitude(double T) const
{
	return get_maxwellian_velocity(T).norm();
}

Vector3d Species::get_maxwellian_velocity(double T) const
{
	double v_th = sqrt(2*K*T/m);

	Vector3d v;
	for(int dim : {X, Y, Z}) {
		v(dim) = sqrt(0.5)*v_th*rng.normal();
	}

	return v;
}

Vector3d Species::get_maxwellian_velocity(const vector<double> T) const
{
	Vector3d v;
	for(int dim : {X, Y, Z}) {
		double v_th = sqrt(2*K*T[dim]/m);
		v(dim) = sqrt(0.5)*v_th*rng.normal();
	}

	return v;
}

void Species::add_particle(const Vector3d &x, const Vector3d &v)
{
	add_particle(x, v, domain.get_time_step(), w_mp0);
}

void Species::add_particle(const Vector3d &x, const Vector3d &v, double dt)
{
	add_particle(x, v, dt, w_mp0);
}

void Species::add_particle(const Vector3d &x, const Vector3d &v, double dt, double w_mp)
{
	Vector3d l = domain.x_to_l(x);
	Vector3d E_p = domain.gather(domain.E, l);
	Vector3d dv = q/m*E_p*0.5*domain.get_time_step();
	particles.push_back(Particle(x, v - dv, dt, w_mp));
}

void Species::add_cold_box(const Vector3d &x1, const Vector3d &x2, double n,
		const Vector3d &v_drift)
{
	double V_box = (x2 - x1).prod();
	double n_real = n*V_box;
	int n_sim = (int)(n_real/w_mp0);

	for (int p = 0; p < n_sim; ++p) {
		Vector3d x;
		do {
			x = x1.array() + rng(3).array()*(x2 - x1).array();
		} while(!domain.is_inside(x));
		add_particle(x, v_drift);
	}
}

void Species::add_warm_box(const Vector3d &x1, const Vector3d &x2, double n,
		const Vector3d &v_drift, const vector<double> T)
{
	double V_box = (x2 - x1).prod();
	double n_real = n*V_box;
	int n_sim = (int)(n_real/w_mp0);

	for (int p = 0; p < n_sim; ++p) {
		Vector3d x;
		do {
			x = x1.array() + rng(3).array()*(x2 - x1).array();
		} while(!domain.is_inside(x));

		Vector3d v_M;
		if (T.size() == 1) {
			v_M = get_maxwellian_velocity(T[0]);
		} else if (T.size() == 3) {
			v_M = get_maxwellian_velocity(T);
		}

		add_particle(x, v_drift + v_M);
	}
}

void Species::push_particles_leapfrog()
{
	for(Particle &p : particles) {
		Vector3d l = domain.x_to_l(p.x);
		Vector3d E_p = domain.gather(domain.E, l);

		p.v += E_p*(p.dt*q/m);

		int n_bounces = 0;

		while (p.dt > 0 && p.w_mp > 0) {
			Vector3d x_old = p.x;
			p.x += p.v*p.dt;

			if (!domain.is_inside(p.x)) {
				domain.apply_boundary_conditions(*this, x_old, p);
				continue;
			}

			p.dt = 0;

			if (++n_bounces > 10) {
				p.w_mp = 0;
			}
		}

		p.dt += domain.get_time_step();
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

	const int &ni = domain.ni;
	const int &nj = domain.nj;
	const int &nk = domain.nk;

	/* synchronize number density, assuming uniform grid spacing */
	if (domain.is_periodic(Xmin)) {
		for (int j = 0; j < nj; ++j) {
			for (int k = 0; k < nk; ++k) {
				n(domain.at(0, j, k)) = 0.5*(n(domain.at(0, j, k))
						+ n(domain.at(ni - 1, j, k)));
				n(domain.at(ni - 1, j, k)) = n(domain.at(0, j, k));
			}
		}
	}

	if (domain.is_periodic(Ymin)) {
		for (int i = 0; i < ni; ++i) {
			for (int k = 0; k < nk; ++k) {
				n(domain.at(i, 0, k)) = 0.5*(n(domain.at(i, 0, k))
						+ n(domain.at(i, nj - 1, k)));
				n(domain.at(i, nj - 1, k)) = n(domain.at(i, 0, k));
			}
		}
	}

	if (domain.is_periodic(Zmin)) {
		for (int i = 0; i < nj; ++i) {
			for (int j = 0; j < nj; ++j) {
				n(domain.at(i, j, 0)) = 0.5*(n(domain.at(i, j, 0))
						+ n(domain.at(i, j, nk - 1)));
				n(domain.at(i, j, nk - 1)) = n(domain.at(i, j, 0));
			}
		}
	}

	n = n.array()/domain.V_node.array();

	/* do the time averaging if mu > 0 */
	n_mean = mu*n_mean + (1 - mu)*n;
}

void Species::sample_moments()
{
	n_sum.setZero();
	nv_sum.setZero();
	nuu_sum.setZero();
	nvv_sum.setZero();
	nww_sum.setZero();

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
		double n_u = n_sum(u);

		if (n_u <= 0) {
			v_stream.row(u).setZero();
			T(u) = 0;
			continue;
		}

		v_stream.row(u) = nv_sum.row(u)/n_u;

		double u_mean = v_stream(u, X);
		double v_mean = v_stream(u, Y);
		double w_mean = v_stream(u, Z);

		double u2_mean = nuu_sum(u)/n_u;
		double v2_mean = nvv_sum(u)/n_u;
		double w2_mean = nww_sum(u)/n_u;

		double uu = u2_mean - u_mean*u_mean;
		double vv = v2_mean - v_mean*v_mean;
		double ww = w2_mean - w_mean*w_mean;

		T(u) = m/(3*K)*(uu + vv + ww);
	}
}

void Species::calc_macroparticle_count()
{
	mp_count.setZero();
	for(const Particle &p : particles) {
		int c = domain.x_to_c(p.x);
		mp_count(c) += 1;
	}
}
