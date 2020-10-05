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
		Options(const int& depth = 8, const int& boundary = 2, const float& samplesPerNode = 1.5, const bool& density = true) :
			depth(depth), boundary(boundary), samplesPerNode(samplesPerNode), density(density) {}
		//Parameters
		int depth;
		//1-BOUNDARY_FREE 2-BOUNDARY_NEUMANN 3-BOUNDARY_DIRICHLET
		int boundary;
		//[1.0-5.0] for point cloud with less noise and [15.0-20.0] if we have a lot of noise
		float samplesPerNode;
		//Generate the value scalar
		bool density;
	};

	bool Compute(AdaptativeSolvers::Mesh<float>& mesh_in_out, const Options& options);
};
