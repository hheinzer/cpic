#include "domain.hpp"
#include "random.hpp"
#include "const.hpp"
#include "species.hpp"

using namespace std;
using namespace Eigen;
using namespace Const;

RandomNumberGenerator rng;

Domain::Domain(string prefix, int ni, int nj, int nk) :
	prefix{prefix}, ni{ni}, nj{nj}, nk{nk}, nn{ni, nj, nk}, n_nodes{ni*nj*nk},
	n_cells{(ni - 1)*(nj - 1)*(nk - 1)}
{
	cout << "┌───────────────────────────────────────────────┐\n"
	     << "│      CPIC ── C++ Particle in Cell Method      │\n"
	     << "│       Written by Heinz Heinrich Heinzer       │\n"
	     << "└───────────────────────────────────────────────┘\n";

	wtime_start = chrono::high_resolution_clock::now();

	V_node = VectorXd::Zero(n_nodes);
	rho    = VectorXd::Zero(n_nodes);
	phi    = VectorXd::Zero(n_nodes);
	E      = MatrixXd::Zero(n_nodes, 3);
	n_e_BR = VectorXd::Zero(n_nodes);
}

void Domain::set_dimensions(const Vector3d &x_min, const Vector3d &x_max)
{
	this->x_min = x_min;
	this->x_max = x_max;
	this->del_x = (x_max - x_min).array()/(nn.cast<double>().array() - 1);
	calc_node_volume();
}

void Domain::set_bc_at(BoundarySide side, BC bc)
{
	if (bc.field_bc_type == FieldBCtype::Neumann) {
		if ((side == Xmin || side == Xmax)) {
			bc.set_delta(del_x(X));
		} else if (side == Ymin || side == Ymax) {
			bc.set_delta(del_x(Y));
		} else if (side == Zmin || side == Zmax) {
			bc.set_delta(del_x(Z));
		}
	}

	this->bc[side] = make_unique<BC>(bc);
}

double Domain::get_wtime() const
{
	auto wtime_now = chrono::high_resolution_clock::now();
	chrono::duration<double> wtime_delta = wtime_now - wtime_start;
	return wtime_delta.count();
}

bool Domain::is_inside(const Vector3d &x) const
{
	return (x_min.array() <     x.array()).all() &&
		   (    x.array() < x_max.array()).all();
}

double Domain::get_potential_energy() const
{
	return 0.5*EPS0*E.rowwise().squaredNorm().transpose()*V_node;
}

Vector3d Domain::x_to_l(const Vector3d &x) const
{
	return (x - x_min).array()/del_x.array();
}

int Domain::x_to_c(const Vector3d &x) const
{
	Vector3i lInt = x_to_l(x).cast<int>();
	return lInt(X) + lInt(Y)*(ni - 1) + lInt(Z)*(ni - 1)*(nj - 1);
}

void Domain::scatter(VectorXd &f, const Vector3d &l, double value)
{
	int i = (int)l(X);
	double di = l(X) - i;

	int j = (int)l(Y);
	double dj = l(Y) - j;

	int k = (int)l(Z);
	double dk = l(Z) - k;

	f(at(i    ,j    ,k    )) += value*(1 - di)*(1 - dj)*(1 - dk);
	f(at(i + 1,j    ,k    )) += value*(    di)*(1 - dj)*(1 - dk);
	f(at(i    ,j + 1,k    )) += value*(1 - di)*(    dj)*(1 - dk);
	f(at(i + 1,j + 1,k    )) += value*(    di)*(    dj)*(1 - dk);
	f(at(i    ,j    ,k + 1)) += value*(1 - di)*(1 - dj)*(    dk);
	f(at(i + 1,j    ,k + 1)) += value*(    di)*(1 - dj)*(    dk);
	f(at(i    ,j + 1,k + 1)) += value*(1 - di)*(    dj)*(    dk);
	f(at(i + 1,j + 1,k + 1)) += value*(    di)*(    dj)*(    dk);
}

void Domain::scatter(MatrixXd &f, const Vector3d &l, const Vector3d &value)
{
	int i = (int)l(X);
	double di = l(X) - i;

	int j = (int)l(Y);
	double dj = l(Y) - j;

	int k = (int)l(Z);
	double dk = l(Z) - k;

	f.row(at(i    ,j    ,k    )) += value*(1 - di)*(1 - dj)*(1 - dk);
	f.row(at(i + 1,j    ,k    )) += value*(    di)*(1 - dj)*(1 - dk);
	f.row(at(i    ,j + 1,k    )) += value*(1 - di)*(    dj)*(1 - dk);
	f.row(at(i + 1,j + 1,k    )) += value*(    di)*(    dj)*(1 - dk);
	f.row(at(i    ,j    ,k + 1)) += value*(1 - di)*(1 - dj)*(    dk);
	f.row(at(i + 1,j    ,k + 1)) += value*(    di)*(1 - dj)*(    dk);
	f.row(at(i    ,j + 1,k + 1)) += value*(1 - di)*(    dj)*(    dk);
	f.row(at(i + 1,j + 1,k + 1)) += value*(    di)*(    dj)*(    dk);
}

Vector3d Domain::gather(const MatrixXd &f, const Vector3d &l) const
{
	int i = (int)l(X);
	double di = l(X) - i;

	int j = (int)l(Y);
	double dj = l(Y) - j;

	int k = (int)l(Z);
	double dk = l(Z) - k;

	return f.row(at(i    ,j    ,k    ))*(1 - di)*(1 - dj)*(1 - dk)
		 + f.row(at(i + 1,j    ,k    ))*(    di)*(1 - dj)*(1 - dk)
		 + f.row(at(i    ,j + 1,k    ))*(1 - di)*(    dj)*(1 - dk)
		 + f.row(at(i + 1,j + 1,k    ))*(    di)*(    dj)*(1 - dk)
		 + f.row(at(i    ,j,    k + 1))*(1 - di)*(1 - dj)*(    dk)
		 + f.row(at(i + 1,j    ,k + 1))*(    di)*(1 - dj)*(    dk)
		 + f.row(at(i    ,j + 1,k + 1))*(1 - di)*(    dj)*(    dk)
		 + f.row(at(i + 1,j + 1,k + 1))*(    di)*(    dj)*(    dk);
}

void Domain::calc_charge_density(std::vector<Species> &species)
{
	rho.setZero();
	for(const Species &sp : species) {
		if (sp.q == 0) continue;
		rho += sp.q*sp.n;
	}
}

void Domain::apply_boundary_conditions(const Species &sp, const Vector3d &x_old,
		Particle &p) const
{
	for(int dim : {X, Y, Z}) {
		if (p.x(dim) < x_min(dim)) {
			int side = 2*dim;
			Vector3d n = Vector3d::Unit(dim);
			eval_particle_BC(sp, side, x_min, x_old, p, dim, n);
		} else if (x_max(dim) < p.x(dim)) {
			int side = 2*dim + 1;
			Vector3d n = -Vector3d::Unit(dim);
			eval_particle_BC(sp, side, x_max, x_old, p, dim, n);
		}
	}
}

void Domain::eval_field_BC(BoundarySide side, VectorXd &b0, std::vector<T> &coeffs,
		int u, int v) const
{
	switch (bc.at(side)->field_bc_type) {
		case FieldBCtype::Dirichlet:
			coeffs.push_back(T(u, u, 1));
			b0(u) = bc.at(side)->get_value();
			break;
		case FieldBCtype::Neumann:
			coeffs.push_back(T(u, u,  1));
			coeffs.push_back(T(u, v, -1));
			b0(u) = bc.at(side)->get_value();
			break;
	}
}

bool Domain::steady_state(std::vector<Species> &species)
{
	if (is_steady_state) return true;

	/* only check every 10 iterations */
	if (iter%10 != 0) return false;

	double n_tot = 0, I_tot = 0, E_tot = 0;
	for(const Species &sp : species) {
		n_tot += sp.get_sim_count();
		I_tot += sp.get_momentum().norm();
		E_tot += sp.get_kinetic_energy();
	}

	if (	abs((n_tot - prev_n_tot)/prev_n_tot) < tol_steady_state &&
			abs((I_tot - prev_I_tot)/prev_I_tot) < tol_steady_state &&
			abs((E_tot - prev_E_tot)/prev_E_tot) < tol_steady_state)  {
		is_steady_state = true;
		cout << "Steady state reached at iteration " << iter << endl;
	}

	prev_n_tot = n_tot;
	prev_I_tot = I_tot;
	prev_E_tot = E_tot;

	return is_steady_state;
}

void Domain::calc_node_volume()
{
	for (int i = 0; i < ni; ++i) {
		for (int j = 0; j < nj; ++j) {
			for (int k = 0; k < nk; ++k) {
				double V = del_x.prod();
				if (i == 0 || i == ni - 1) V /= 2;
				if (j == 0 || j == nj - 1) V /= 2;
				if (k == 0 || k == nk - 1) V /= 2;
				V_node(at(i, j, k)) = V;
			}
		}
	}
}

void Domain::eval_particle_BC(const Species &sp, int side, const Vector3d &X,
		const Vector3d &x_old, Particle &p, int dim, const Vector3d &n) const
{
	switch (bc.at(side)->particle_bc_type) {
		case ParticleBCtype::Symmetric:
		case ParticleBCtype::Specular:
			p.x(dim) = 2*X(dim) - p.x(dim);
			p.v(dim) *= -1;
			break;
		case ParticleBCtype::Open:
			p.w_mp = 0;
			break;
		case ParticleBCtype::Diffuse: {
				double t = (X(dim) - x_old(dim))/(p.x(dim) - x_old(dim));
				double dt_rem = (1 - t)*p.dt;
				p.dt -= dt_rem;

				p.x = x_old + 0.999*t*(p.x - x_old);
				double v_mag1 = p.v.norm();

				double v_th = sp.get_maxwellian_velocity_magnitude(bc.at(side)->T);
				double v_mag2 = v_mag1 + bc.at(side)->a_th*(v_th - v_mag1);
				p.v = v_mag2*get_diffuse_vector(n);
				break;
			}
	}
}

Vector3d Domain::get_diffuse_vector(const Vector3d &n) const
{
	/* random vector that follows the cosine law */
	double sin_theta = rng();
	double cos_theta = sqrt(1 - sin_theta*sin_theta);
	double psi = 2*PI*rng();

	Vector3d t1;
	if (n.cross(Vector3d::UnitX()).norm() != 0) {
		t1 = n.cross(Vector3d::UnitX());
	} else {
		t1 = n.cross(Vector3d::UnitY());
	}
	Vector3d t2 = n.cross(t1);

	return sin_theta*(cos(psi)*t1 + sin(psi)*t2) + cos_theta*n;
}

void Domain::check_formulation(double n_e, double T_e, const std::vector<double> &z_i,
		const std::vector<double> &n_i, const std::vector<double> &T_i) const
{
	cout << "Formulation Check:" << endl;

	double dx = del_x.maxCoeff();
	cout << "  Grid spacing:                   "
		<< dx << " m" << endl;
	cout << "  Timestep:                       "
		<< dt << " s" << endl;

	double n_T_i_sum = 0;
	for(size_t i = 0; i < z_i.size(); ++i)
		n_T_i_sum += z_i[i]*z_i[i]*n_i[i]/T_i[i];
	double lambda_D = sqrt(EPS0*K/(QE*QE)/(n_e/T_e + n_T_i_sum));
	cout << "  Debye length:                   "
		<< lambda_D << " m" << endl;

	double omega_p = sqrt(n_e*QE*QE/(EPS0*ME));
	cout << "  Plasma frequency:               "
		<< omega_p << " rad/s" << endl;

	cout << "  Spatial stability:  Δx < λ_D:   "
		 << (dx < lambda_D ? "true" : "false") << endl;

	cout << "  Temporal stability: Δt < 1/ω_p: "
		 << (dt < 1/omega_p ? "true" : "false") << endl << endl;
}

void Domain::print_info(std::vector<Species> &species) const
{
	cout << "iter:" << setw(6) << iter;

	double v_max = 0;
	for(const Species &sp : species) {
		cout << "  " << sp.name << ":" << setw(6) << sp.get_sim_count();

		#pragma omp parallel for
		for(const Particle &p : sp.particles)
			if (p.v.norm() > v_max)
				v_max = p.v.norm();
	}

	cout << "  CFL: " << setprecision(3) << v_max*dt/del_x.minCoeff() << endl;
}

void Domain::write_statistics(std::vector<Species> &species)
{
	if (!stats.is_open()) {
		stats.open(prefix + "_statistics.csv");
		stats << "iter,time,wtime";
		for(const Species &sp : species) {
			stats << ",n_sim." << sp.name
				  << ",n_real." << sp.name
				  << ",Ix." << sp.name
				  << ",Iy." << sp.name
				  << ",Iz." << sp.name
				  << ",E_kin." << sp.name;
		}
		stats << ",E_pot,E_tot" << endl;
	}

	stats << iter << "," << time << "," << get_wtime() << ",";

	double E_tot = 0;
	for(const Species &sp : species) {
		double E_kin = sp.get_kinetic_energy();
		Vector3d I = sp.get_momentum();
		stats << sp.get_sim_count() << ","
			  << sp.get_real_count() << ","
			  << I(X) << ","
			  << I(Y) << ","
			  << I(Z) << ","
			  << E_kin << ",";

		E_tot += E_kin;
	}

	double E_pot = get_potential_energy();
	stats << E_pot << "," << E_tot + E_pot << "\n";

	if (iter%25 == 0) stats.flush();
}

void Domain::save_fields(std::vector<Species> &species) const
{
	stringstream ss;
	ss << prefix << "_" << setfill('0') << setw(6) << get_iter() << ".vti";

	ofstream out(ss.str());
	if (!out.is_open()) {
		cerr << "Could not open '" << ss.str() << "'" << endl;
		exit(EXIT_FAILURE);
	}

	out << "<VTKFile type=\"ImageData\">\n";
	out << "<ImageData Origin=\"" << x_min.transpose() << "\" ";
	out << "Spacing=\"" << del_x.transpose() << "\" ";
	out << "WholeExtent=\"0 " << ni - 1
					 << " 0 " << nj - 1
					 << " 0 " << nk - 1 << "\">\n";

	out << "<PointData>\n";

	out << "<DataArray Name=\"V_node\" NumberOfComponents=\"1\" "
		<< "format=\"ascii\" type=\"Float64\">\n";
	out << V_node;
	out << "</DataArray>\n";

	out << "<DataArray Name=\"rho\" NumberOfComponents=\"1\" "
		<< "format=\"ascii\" type=\"Float64\">\n";
	out << rho;
	out << "</DataArray>\n";

	out << "<DataArray Name=\"phi\" NumberOfComponents=\"1\" "
		<< "format=\"ascii\" type=\"Float64\">\n";
	out << phi;
	out << "</DataArray>\n";

	out << "<DataArray Name=\"E\" NumberOfComponents=\"3\" "
		<< "format=\"ascii\" type=\"Float64\">\n";
	out << E;
	out << "</DataArray>\n";

	out << "<DataArray Name=\"n.fluid_e-\" NumberOfComponents=\"1\" "
		<< "format=\"ascii\" type=\"Float64\">\n";
	out << n_e_BR;
	out << "</DataArray>\n";

	for (const Species &sp : species) {
		out << "<DataArray Name=\"n." << sp.name
			<< "\" NumberOfComponents=\"1\" format=\"ascii\" "
			<< "type=\"Float64\">\n";
		out << sp.n;
		out << "</DataArray>\n";

		out << "<DataArray Name=\"n_mean." << sp.name
			<< "\" NumberOfComponents=\"1\" format=\"ascii\" "
			<< "type=\"Float64\">\n";
		out << sp.n_mean;
		out << "</DataArray>\n";

		out << "<DataArray Name=\"v_stream." << sp.name
			<< "\" NumberOfComponents=\"3\" format=\"ascii\" "
			<< "type=\"Float64\">\n";
		out << sp.v_stream;
		out << "</DataArray>\n";

		out << "<DataArray Name=\"T." << sp.name
			<< "\" NumberOfComponents=\"1\" format=\"ascii\" "
			<< "type=\"Float64\">\n";
		out << sp.T;
		out << "</DataArray>\n";
	}

	out << "</PointData>\n";

	out << "<CellData>\n";
	for (const Species &sp : species) {
		out << "<DataArray Name=\"mp_count." << sp.name
			<< "\" NumberOfComponents=\"1\" format=\"ascii\" "
			<< "type=\"Float64\">\n";
		out << sp.mp_count;
		out << "</DataArray>\n";
	}
	out << "</CellData>\n";

	out << "</ImageData>\n";
	out << "</VTKFile>\n";

	out.close();
}

void Domain::save_particles(std::vector<Species> &species, int n_particles) const
{
	for(const Species &sp : species) {
		stringstream ss;
		ss << prefix << "_" << sp.name << "_particles"
			<< "_" << setfill('0') << setw(6) << get_iter() << ".vtp";

		ofstream out(ss.str());
		if (!out.is_open()) {
			cerr << "Could not open '" << ss.str() << "'" << endl;
			exit(EXIT_FAILURE);
		}

		double dp = n_particles/(double)sp.get_sim_count();
		double np = 0;
		vector<const Particle *> p_out;
		for(const Particle &p : sp.particles) {
			np += dp;
			if (np > 1) {
				p_out.emplace_back(&p);
				np -= 1;
			}
		}

		out << "<?xml version=\"1.0\"?>\n";
		out << "<VTKFile type=\"PolyData\" version=\"0.1\" "
			<< "byte_order=\"LittleEndian\">\n";
		out << "<PolyData>\n";
		out << "<Piece NumberOfPoints=\"" << p_out.size()
			<< "\" NumberOfVerts=\"0\" NumberOfLines=\"0\" ";
		out << "NumberOfStrips=\"0\" NumberOfCells=\"0\">\n";

		out << "<Points>\n";
		out << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" "
			<< "format=\"ascii\">\n";
		for (const Particle *p : p_out)
			out << p->x.transpose() << "\n";
		out << "</DataArray>\n";
		out << "</Points>\n";

		out << "<PointData>\n";
		out << "<DataArray Name=\"v." << sp.name
			<< "\" type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
		for (const Particle *p : p_out)
			out << p->v.transpose() << "\n";
		out << "</DataArray>\n";
		out << "</PointData>\n";

		out << "</Piece>\n";
		out << "</PolyData>\n";
		out << "</VTKFile>\n";

		out.close();
	}
}

void Domain::save_velocity_histogram(std::vector<Species> &species) const
{
	for(const Species &sp :species) {
		stringstream ss;
		ss << prefix << "_v." << sp.name << "_histogram"
			<< "_" << setfill('0') << setw(6) << get_iter() << ".csv";

		ofstream out(ss.str());
		if (!out.is_open()) {
			cerr << "Could not open '" << ss.str() << "'" << endl;
			exit(EXIT_FAILURE);
		}

		const int n_bins = 100;

		Vector3d v_min = Vector3d::Zero();
		Vector3d v_max = Vector3d::Zero();
		double v_mag_max = 0;

		if (sp.particles.size() > 0) {
			v_min = sp.particles[0].v;
			v_max = sp.particles[0].v;
			v_mag_max = sp.particles[0].v.norm();
		}

		for(const Particle &p : sp.particles) {
			v_min = v_min.cwiseMin(p.v);
			v_max = v_max.cwiseMax(p.v);
			v_mag_max = max(v_mag_max, p.v.norm());
		}

		const double eps = 0.01;
		v_min.array() -= eps;
		v_max.array() += eps;
		v_mag_max += eps;

		Vector4d dv;
		dv << (v_max - v_min)/(n_bins - 1), v_mag_max/(n_bins - 1);

		MatrixXd bins = MatrixXd::Zero(n_bins, 4);
		for(const Particle &p : sp.particles) {
			for(int dim : {X, Y, Z}) {
				bins(round((p.v(dim) - v_min(dim))/dv(dim)), dim) += 1;
			}
			bins(round(p.v.norm()/dv(W)), W) += 1;
		}

		bins /= sp.particles.size();

		out << "v_x,f(v_x),v_y,f(v_y),v_z,f(v_z),v_mag,f(v_mag)\n";
		for (int b = 0; b < n_bins; ++b) {
			for(int dim : {X, Y, Z})
				out << v_min(dim) + b*dv(dim) << "," << bins(b, dim) << ",";
			out << b*dv(W) << "," << bins(b, W) << "\n";
		}

		out.close();
	}
}
