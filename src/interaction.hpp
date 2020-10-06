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

class DSMCneutral : public Interaction {
	public:
		using Vector3d = Eigen::Vector3d;
		using RefPropMap = std::map<std::string, double>;

		DSMCneutral(Domain &domain, Species &species);

		void apply(double dt) override;

	private:
		Domain &domain;
		Species &species;

		const RefPropMap d_ref_map = {	/* [m] reference diameter */
			{"O", 4.07e-10}
		};
		const RefPropMap T_ref_map = {	/* [K] reference temperature */
			{"O", 273.15}
		};
		const RefPropMap omega_map = {	/* [-] viscosity index */
			{"O", 0.77}
		};

		double d_ref, T_ref, omega;

		int n_cells;

		double V;		/* [m^3] cell volume */
		double w_mp;	/* [-] macroparticle weight */

		double sigma_vr_max = 1e-14;

		double mr, m;	/* [kg] reduced mass, species mass */
		double df;		/* [m^2] VHS diameter factor squared */

		double sigma(double v_r) const;

		void collide(Vector3d &v1, Vector3d &v2, double m1, double m2) const;
};

#endif
