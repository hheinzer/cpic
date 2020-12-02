#ifndef RANDOM_HPP
#define RANDOM_HPP

#include <random>
#include <Eigen/Dense>

class RandomNumberGenerator {
	public:
		using VectorXd = Eigen::VectorXd;
		using MatrixXd = Eigen::MatrixXd;

		RandomNumberGenerator() :
			gen{std::random_device()()}, uni{0.0, 1.0}, nrm{0.0, 1.0} {}

		double operator()() {return uni(gen);}

		VectorXd operator()(int ni) {
			return VectorXd::NullaryExpr(ni, [&](){return uni(gen);});
		}

		MatrixXd operator()(int ni, int nj) {
			return MatrixXd::NullaryExpr(ni, nj, [&](){return uni(gen);});
		}

		double normal() {return nrm(gen);}

		std::mt19937& get_gen() {return gen;}

	private:
		std::mt19937 gen;
		std::uniform_real_distribution<double> uni;
		std::normal_distribution<double> nrm;
};

extern RandomNumberGenerator rng;

#endif
