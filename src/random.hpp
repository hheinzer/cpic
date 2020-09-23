#ifndef RANDOM_HPP
#define RANDOM_HPP

#include <random>
#include <Eigen/Dense>

class RandomNumberGenerator {
	public:
		using VectorXd = Eigen::VectorXd;
		using MatrixXd = Eigen::MatrixXd;

		RandomNumberGenerator() : gen{std::random_device()()}, dist{0, 1} {}

		double operator()() {return dist(gen);}

		VectorXd operator()(int ni) {
			return VectorXd::NullaryExpr(ni, [&](){return dist(gen);});
		}

		MatrixXd operator()(int ni, int nj) {
			return MatrixXd::NullaryExpr(ni, nj, [&](){return dist(gen);});
		}

	private:
		std::mt19937 gen;
		std::uniform_real_distribution<double> dist;
};

extern RandomNumberGenerator rng;

#endif
