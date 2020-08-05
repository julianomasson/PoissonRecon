#ifdef SSDRECON_EXPORTS
#define SSDRECON_API __declspec(dllexport)
#else
#define SSDRECON_API __declspec(dllimport)
#endif

#pragma once

#include <vector>

namespace SSD
{
	template <typename T>
	struct Point
	{
		T xyz[3];
		T normal[3];
		unsigned char color[3];
		T value;
	};

	struct Face
	{
		unsigned int point_indices[3];
	};

	template <typename T>
	struct Mesh
	{
		std::vector< Point<T> > points;
		std::vector< Face > faces;
		//Parameters
		bool has_normal = false;
		bool has_color = false;
	};

	// Define the poisson options
	struct Options
	{
		Options() {}
		Options(const int& depth = 8, const int& boundary = 2, const bool& density = true) : depth(depth), boundary(boundary), density(density) {}
		//Parameters
		int depth;
		//1-BOUNDARY_FREE 2-BOUNDARY_NEUMANN 3-BOUNDARY_DIRICHLET
		int boundary;
		//Generate the value scalar
		bool density;
	};
}

class SSDRECON_API SSDRecon
{
public:
	SSDRecon(void) {};
	template <typename T>
	bool compute(SSD::Mesh<T>& mesh_in_out, const SSD::Options& options);
};

