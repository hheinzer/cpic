#ifndef OBJECT_HPP
#define OBJECT_HPP

#include <map>
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include "domain.hpp"

enum class ElementType {Triangle = 2};

struct Node {
	using Vector3d = Eigen::Vector3d;

	Node(double x, double y, double z) : pos{x, y, z} {}

	const Vector3d pos;
};

struct Triangle {
	using Vector3i = Eigen::Vector3i;
	using Vector3d = Eigen::Vector3d;

	Triangle(const Vector3i &con, const Vector3d &x_c, const Vector3d &n) :
		con{con}, x_c{x_c}, n{n} {}

	const Vector3i con;
	const Vector3d x_c;
	const Vector3d n;
};

class SurfaceMesh {
	public:
		using Vector3d = Eigen::Vector3d;

		SurfaceMesh(std::string file_name) {load_msh(file_name);}

		std::vector<Node> nodes;
		std::vector<Triangle> trias;

	private:
		std::map<int, int> lut_nodes;
		std::map<int, int> lut_triangles;

		void load_msh(std::string file_name);

		void load_msh_nodes(std::ifstream &msh);

		void load_msh_elements(std::ifstream &msh);
};

class Object {
	public:
		using Vector3d = Eigen::Vector3d;
		using Intersection = std::pair<double, Vector3d>;

		Object(Domain &domain, std::string mesh_file_name, const Vector3d &x_o,
				double phi, double T);

		bool is_inside(const Vector3d &x) const;

		bool is_inside(int u) const;

		Intersection get_intersection(const Vector3d &x1, const Vector3d &x2) const;

		Intersection get_diffuse_intersection(const Vector3d &x1, const Vector3d &x2) const;

		const Vector3d x_o;	/* [m] origin offset vector */
		const double phi;	/* [V] surface potential */
		const double T;		/* [K] surface temperature */

	private:
		Domain &domain;
		SurfaceMesh mesh;
};

#endif
