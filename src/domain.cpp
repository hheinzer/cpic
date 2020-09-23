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
	rho = VectorXd::Zero(n_nodes);
	phi = VectorXd::Zero(n_nodes);
	E = MatrixXd::Zero(n_nodes, 3);
}

void Domain::set_dimensions(const Vector3d &x_min, const Vector3d &x_max)
{
	this->x_min = x_min;
	this->x_max = x_max;
	this->del_x = (x_max - x_min).array()/(nn.cast<double>().array() - 1);
	calc_node_volume();
}

void Domain::set_boundary_condition(BoundarySide side, BC bc)
{
	if (side == Xmin || side == Xmax) {
		bc.set_delta(del_x(X));
	} else if (side == Ymin || side == Ymax) {
		bc.set_delta(del_x(Y));
	} else if (side == Zmin || side == Zmax) {
		bc.set_delta(del_x(Z));
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

void Domain::scatter(MatrixXd &f, const Vector3d &l, const VectorXd &value)
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
		if (sp.rho_s == 0) continue;
		rho += sp.rho_s*sp.n;
	}
}

void Domain::apply_boundary_conditions(Particle &p)
{
	for(int dim : {X, Y, Z}) {
		if (p.x(dim) < x_min(dim)) {
			int side = 2*dim;
			ParticleBCtype type = bc.at(side)->particle_bc_type;
			eval_particle_BC(type, x_min(dim), p.x(dim), p.v(dim), p.w_mp);
		} else if (x_max(dim) < p.x(dim)) {
			int side = 2*dim + 1;
			ParticleBCtype type = bc.at(side)->particle_bc_type;
			eval_particle_BC(type, x_max(dim), p.x(dim), p.v(dim), p.w_mp);
		}
	}
}

void Domain::eval_field_BC(BoundarySide side, VectorXd &b0,
		std::vector<T> &coeffs, int u, int v)
{
	switch (bc.at(side)->field_bc_type) {
		case Dirichlet:
			coeffs.push_back(T(u, u, 1));
			b0(u) = bc.at(side)->get_value();
			break;
		case Neumann:
			coeffs.push_back(T(u, u,  1));
			coeffs.push_back(T(u, v, -1));
			b0(u) = bc.at(side)->get_value();
			break;
	}
}

bool Domain::steady_state(std::vector<Species> &species)
{
	if (is_steady_state) return true;

	double m_tot = 0, I_tot = 0, E_tot = 0;
	for(const Species &sp : species) {
		m_tot += sp.get_real_count();
		I_tot += sp.get_momentum().norm();
		E_tot += sp.get_kinetic_energy();
	}

	if (	abs((m_tot - prev_m_tot)/prev_m_tot) < tol_steady_state &&
			abs((I_tot - prev_I_tot)/prev_I_tot) < tol_steady_state &&
			abs((E_tot - prev_E_tot)/prev_E_tot) < tol_steady_state)  {
		is_steady_state = true;
		cout << "Steady state reached at iteration " << iter << endl;
	}

	prev_m_tot = m_tot;
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

void Domain::eval_particle_BC(ParticleBCtype type, const double &X, double &x,
		double &v, double &w_mp)
{
	switch (type) {
		case Reflective:
			x = 2*X - x;
			v *= -1;
			break;
		case Open:
			w_mp = 0;
			break;
	}
}

void Domain::print_info(std::vector<Species> &species) const
{
	cout << "iter: " << setw(6) << iter;
	for(const Species &sp : species)
		cout << "\t" << sp.name << ": " << setw(6) << sp.get_sim_count();
	cout << endl;
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

	for (Species &sp : species) {
		sp.calc_gas_properties();

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
	for (Species &sp : species) {
		sp.calc_macroparticle_density();

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

	if (!steady_state()) {
		for(Species &sp : species)
			sp.clear_samples();
	}
}
