#include "source.hpp"

using namespace std;
using namespace Eigen;

ColdBeam::ColdBeam(Species &species, Domain &domain, const Vector3d &x1,
		const Vector3d &x2, const Vector3d &v0, double n) :
	species{species}, domain{domain}, x1{x1}, x2{x2}, v0{v0}, n{n}
{
	dx = x2 - x1;
	if (dx(X) == 0) {
		fix_dir        = X;
		random_dirs[0] = Y;
		random_dirs[1] = Z;
		A = dx(Y)*dx(Z);
	} else if (dx(Y) == 0) {
		random_dirs[0] = X;
		fix_dir        = Y;
		random_dirs[1] = Z;
		A = dx(X)*dx(Z);
	} else if (dx(Y) == 0) {
		random_dirs[0] = X;
		random_dirs[1] = Y;
		fix_dir        = Z;
		A = dx(X)*dx(Y);
	} else {
		cerr << "Source definition error: x1 and x2 don't define a plane" << endl;
		exit(EXIT_FAILURE);
	}
}

void ColdBeam::sample()
{
	double n_real = n*v0.norm()*A*domain.get_time_step();
	int n_sim = (int)(n_real/species.w_mp0 + rng());

	MatrixXd new_x(n_sim, 3), new_v(n_sim, 3);

	for (int &dir : random_dirs)
		new_x.col(dir) = x1(dir) + rng(n_sim).array()*dx(dir);
	new_x.col(fix_dir).setConstant(x1(fix_dir));

	for (int dir : {X, Y, Z})
		new_v.col(dir).setConstant(v0(dir));

	for (int p = 0; p < n_sim; ++p) {
		species.add_particle(new_x.row(p), new_v.row(p));
	}
}
