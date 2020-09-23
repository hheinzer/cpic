#ifndef CONST_HPP
#define CONST_HPP

namespace Const
{
	constexpr double EPS0  = 8.85418782e-12;	/* [C/(V*m)] vacuum permittivity */
	constexpr double QE	   = 1.602176565e-19;	/* [C] electron charge */
	constexpr double AMU   = 1.660538921e-27;	/* [kg] atomic mass unit */
	constexpr double ME	   = 9.10938215e-31;	/* [kg] electron mass */
	constexpr double K	   = 1.380648e-23;		/* [J/K] Boltzmann constant */
	constexpr double PI	   = 3.141592653589793; /* [rad] pi */
	constexpr double EvToK = QE/K;				/* 1eV in K ~ 11604 */
}

#endif
