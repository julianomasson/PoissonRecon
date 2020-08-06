#ifdef SSDRECON_EXPORTS
#define SSDRECON_API __declspec(dllexport)
#else
#define SSDRECON_API __declspec(dllimport)
#endif
#pragma once

#include <vector>
#include "MeshData.h"

class SSDRECON_API SSDRecon
{
public:
	SSDRecon(void) {};

	// Define the SSD options
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

	template <typename T>
	bool compute(AdaptativeSolvers::Mesh<T>& mesh_in_out, const Options& options);
};

