#include "solver.hpp"
#include "const.hpp"
#include "domain.hpp"

using namespace std;
using namespace Eigen;
using namespace Const;

Solver::Solver(Domain &domain, int iter_max, double tol) :
	domain{domain}, iter_max{iter_max}, tol{tol}
{
	Vector3d del_x = domain.get_del_x();
	Vector3d del_x_2q = 1.0/del_x.array().pow(2);

	n_nodes = domain.n_nodes;
	vector<T> coeffs;

	b0 = VectorXd::Zero(n_nodes);
	is_regular = VectorXd::Zero(n_nodes);

	const int &ni = domain.ni;
	const int &nj = domain.nj;
	const int &nk = domain.nk;

	for (int i = 0; i < ni; ++i) {
		for (int j = 0; j < nj; ++j) {
			for (int k = 0; k < nk; ++k) {
				int u = at(i,j,k);

				if (i == 0) {
					domain.eval_field_BC(Xmin, b0, coeffs, u, at(i + 1,j,k));
				} else if (i == domain.ni - 1) {
					domain.eval_field_BC(Xmax, b0, coeffs, u, at(i - 1,j,k));
				} else if (j == 0) {
					domain.eval_field_BC(Ymin, b0, coeffs, u, at(i,j + 1,k));
				} else if (j == domain.nj - 1) {
					domain.eval_field_BC(Ymax, b0, coeffs, u, at(i,j - 1,k));
				} else if (k == 0) {
					domain.eval_field_BC(Zmin, b0, coeffs, u, at(i,j,k + 1));
				} else if (k == domain.nk - 1) {
					domain.eval_field_BC(Zmax, b0, coeffs, u, at(i,j,k - 1));
				} else {
					coeffs.push_back(T(u, u,	  -2*(del_x_2q.sum())));

					coeffs.push_back(T(u, at(i + 1,j,k), del_x_2q(X)));
					coeffs.push_back(T(u, at(i - 1,j,k), del_x_2q(X)));

					coeffs.push_back(T(u, at(i,j + 1,k), del_x_2q(Y)));
					coeffs.push_back(T(u, at(i,j - 1,k), del_x_2q(Y)));

					coeffs.push_back(T(u, at(i,j,k + 1), del_x_2q(Z)));
					coeffs.push_back(T(u, at(i,j,k - 1), del_x_2q(Z)));

					is_regular(u) = 1;
				}
			}
		}
	}

	A.resize(n_nodes, n_nodes);
	A.setFromTriplets(coeffs.begin(), coeffs.end());

	solver.setMaxIterations(iter_max);
	solver.setTolerance(tol);
	solver.compute(A);

	if (solver.info() != Success) {
		cerr << "Solver failed to decompose Matrix!" << endl;
		exit(EXIT_FAILURE);
	}
}

void Solver::set_reference_values(double phi0, double Te0, double n0)
{
	this->phi0 = phi0;
	this->Te0 = Te0;
	this->n0 = n0;
}

void Solver::calc_potential()
{
	VectorXd &rho = domain.rho;
	VectorXd &phi = domain.phi;

	VectorXd b = b0.array() - (rho/EPS0).array()*is_regular.array();

	domain.phi = solver.solveWithGuess(b, phi);

	if (solver.info() != Success) {
		cerr << "Solver failed to find a solution!" << endl;
		exit(EXIT_FAILURE);
	}
}

void Solver::calc_potential_non_linear()
{
	VectorXd &rho = domain.rho;
	VectorXd &phi = domain.phi;
	VectorXd &n_e = domain.n_e;

	VectorXd b = b0.array() - (rho/EPS0).array()*is_regular.array();

	VectorXd del_phi = VectorXd::Zero(n_nodes);

	for (int iter = 0; iter < newton_iter_max; ++iter) {
		VectorXd R = A*phi - b;

		R.array() -= (QE/EPS0*n0*exp((phi.array() - phi0)/Te0))
			*is_regular.array();

		SpMat J = A;

		J.diagonal().array() -= (QE*n0/(EPS0*Te0)*exp((phi.array() - phi0)/Te0))
			*is_regular.array();

		del_phi = solver.factorize(J).solveWithGuess(R, del_phi);

		if (solver.info() != Success) {
			cerr << "Solver failed to find a solution!" << endl;
			exit(EXIT_FAILURE);
		}

		phi -= del_phi;

		if (del_phi.norm() < newton_tol) {
			n_e = n0 * exp((phi.array() - phi0)/Te0);
			return;
		}
	}

	cerr << "Newton Sover failed to converge!" << endl;
	exit(EXIT_FAILURE);
}

void Solver::calc_electric_field()
{
	const int &ni = domain.ni;
	const int &nj = domain.nj;
	const int &nk = domain.nk;

	VectorXd &phi = domain.phi;
	MatrixXd &E  = domain.E;

	Vector3d del_x = domain.get_del_x();
	double dx2 = 2*del_x(X);
	double dy2 = 2*del_x(Y);
	double dz2 = 2*del_x(Z);

	for (int i = 0; i < ni; ++i) {
		for (int j = 0; j < nj; ++j) {
			for (int k = 0; k < nk; ++k) {
				int u = at(i, j, k);

				if (i == 0) {
					E(u,X) = -(-3*phi(at(i,j,k)) + 4*phi(at(i + 1,j,k)) - phi(at(i + 2,j,k)))/dx2;
				} else if (i == ni - 1) {
					E(u,X) = -(phi(at(i - 2,j,k)) - 4*phi(at(i - 1,j,k)) + 3*phi(at(i,j,k)))/dx2;
				} else {
					E(u,X) = -(phi(at(i + 1,j,k)) - phi(at(i - 1,j,k)))/dx2;
				}

				if (j == 0) {
					E(u,Y) = -(-3*phi(at(i,j,k)) + 4*phi(at(i,j + 1,k)) - phi(at(i,j + 2,k)))/dy2;
				} else if (j == nj - 1) {
					E(u,Y) = -(phi(at(i,j - 2,k)) - 4*phi(at(i,j - 1,k)) + 3*phi(at(i,j,k)))/dy2;
				} else {
					E(u,Y) = -(phi(at(i,j + 1,k)) - phi(at(i,j - 1,k)))/dy2;
				}

				if (k == 0) {
					E(u,Z) = -(-3*phi(at(i,j,k)) + 4*phi(at(i,j,k + 1)) - phi(at(i,j,k + 2)))/dz2;
				} else if (k == nk - 1) {
					E(u,Z) = -(phi(at(i,j,k - 2)) - 4*phi(at(i,j,k - 1)) + 3*phi(at(i,j,k)))/dz2;
				} else {
					E(u,Z) = -(phi(at(i,j,k + 1)) - phi(at(i,j,k - 1)))/dz2;
				}
			}
		}
	}
}
