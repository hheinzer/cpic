#include "object.hpp"

using namespace std;
using namespace Eigen;

SurfaceMesh::SurfaceMesh(std::string file_name, const Vector3d &x_o) :
	x_o{x_o}
{
	load_msh(file_name);
}

void SurfaceMesh::load_msh(string file_name)
{
	ifstream msh(file_name);
	if (!msh) {
		cerr << "Could not open '" << file_name << "'" << endl;
		exit(EXIT_FAILURE);
	}

	string line;
	stringstream ss;

	while (getline(msh, line)) {
		if (line != "$MeshFormat") continue;
		getline(msh, line); ss << line;
		float msh_version;
		ss >> msh_version;
		if (msh_version < 4.0) {
			cerr << "Wrong MeshFormat! Only 4.0 and up supported!" << endl;
			exit(EXIT_FAILURE);
		} else {
			cout << "Reading surface mesh file '" << file_name << "':" << endl;
		}
		break;
	}

	while (getline(msh, line)) {
		if (line == "$Nodes") {
			load_msh_nodes(msh);
		} else if (line == "$Elements") {
			load_msh_elements(msh);
		}
	}
}

void SurfaceMesh::load_msh_nodes(ifstream &msh)
{
	string line;
	getline(msh, line);
	stringstream ss(line);

	int num_entity_blocks, num_nodes;
	ss >> num_entity_blocks >> num_nodes;

	float skip;

	for (int i = 0; i < num_entity_blocks; ++i) {
		getline(msh, line);
		stringstream ss(line);
		int num_nodes_in_block;
		ss >> skip >> skip >> skip >> num_nodes_in_block;

		vector<int> id;
		for (int j = 0; j < num_nodes_in_block; ++j) {
			getline(msh, line);
			stringstream ss(line);
			int nodeTag;
			ss >> nodeTag;
			id.emplace_back(nodeTag);
		}

		for (int j = 0; j < num_nodes_in_block; ++j) {
			getline(msh, line);
			stringstream ss(line);
			double x, y, z;
			ss >> x >> y >> z;
			lut_nodes[id[j]] = nodes.size();
			nodes.emplace_back(Node(x + x_o(X), y + x_o(Y), z + x_o(Z)));
		}
	}

	if (num_nodes == (int)nodes.size()) {
		cout << "Number of nodes: " << num_nodes << endl;
	} else {
		cerr << "Error while reading nodes!" << endl;
		exit(EXIT_FAILURE);
	}
}

void SurfaceMesh::load_msh_elements(ifstream &msh)
{
	string line;
	getline(msh, line);
	stringstream ss(line);

	int num_entity_blocks, num_elements;
	ss >> num_entity_blocks >> num_elements;

	float skip;

	for (int i = 0; i < num_entity_blocks; ++i) {
		getline(msh, line);
		stringstream ss(line);
		int element_type, num_elements_in_block;
		ss >> skip >> skip >> element_type >> num_elements_in_block;

		switch ((ElementType)element_type) {
			case ElementType::Triangle:
				for (int j = 0; j < num_elements_in_block; ++j) {
					getline(msh, line);
					stringstream ss(line);

					int id, n1, n2, n3;
					ss >> id >> n1 >> n2 >> n3;
					lut_triangles[id] = trias.size();

					Vector3i con(lut_nodes[n1], lut_nodes[n2], lut_nodes[n3]);

					Vector3d v0 = nodes[con(X)].pos;
					Vector3d v1 = nodes[con(Y)].pos;
					Vector3d v2 = nodes[con(Z)].pos;

					Vector3d x_c = (v0 + v1 + v2)/3;

					Vector3d n = (v1 - v0).cross(v2 - v0).normalized();

					trias.emplace_back(Triangle(con, x_c, n));
				}
				break;
		}
	}

	cout << "Number of triangles read: " << trias.size() << endl;
}

Object::Object(Domain &domain, std::string mesh_file_name, const Vector3d &x_o,
		double phi, double T) :
	phi{phi}, T{T}, mesh{mesh_file_name, x_o}, domain{domain}
{

}

bool Object::is_inside(const Vector3d &x) const
{

}

bool Object::is_inside(int u) const
{

}

Object::Intersection Object::get_intersection(const Vector3d &x1,
		const Vector3d &x2) const
{

}

Object::Intersection Object::get_diffuse_intersection(const Vector3d &x1,
		const Vector3d &x2) const
{

}
