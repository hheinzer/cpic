#ifndef INTERACTION_HPP
#define INTERACTION_HPP

#include <map>
#include <string>
#include <Eigen/Dense>
#include "domain.hpp"
#include "species.hpp"

class Interaction {
	public:
		virtual void apply(double dt) = 0;

		virtual ~Interaction() {}
};

class DSMC : public Interaction {
	public:
		using Vector3d = Eigen::Vector3d;
		using RefPropMap = std::map<std::string, double>;

		DSMC(Domain &domain, Species &sp1, Species &sp2);

		void apply(double dt) override;

	private:
		Domain &domain;
		Species &sp1, &sp2;

		const RefPropMap d_ref = {	/* [m] reference diameter */
			{"O", 4.07e-10}
		};
		const RefPropMap T_ref = {	/* [K] reference temperature */
			{"O", 273.15}
		};
		const RefPropMap omega = {	/* [-] viscosity index */
			{"O", 0.77}
		};

		double d_ref1, d_ref2, T_ref1, T_ref2, omega1, omega2;

		int n_cells;

		double V;		/* [m^3] cell volume */
		double w_mp;	/* [-] macroparticle weight */

		double sigma_vr_max = 1e-14;

		double mr, m1, m2;	/* [kg] reduced mass, species mass */
		double df1, df2;	/* [m^2] VHS diameter factor squared */

		double sigma(double v_r) const;

		void collide(Vector3d &v1, Vector3d &v2, double m1, double m2) const;
};

#endif
