#ifdef POISSONRECON_EXPORTS
#define POISSONRECON_API __declspec(dllexport)
#else
#define POISSONRECON_API __declspec(dllimport)
#endif
#pragma once

#include <vector>
#include "MeshData.h"

class POISSONRECON_API PoissonRecon
{
public:
	PoissonRecon(void) {};

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

	template <typename T>
	bool compute(AdaptativeSolvers::Mesh<T>& mesh_in_out, const Options& options);
};
