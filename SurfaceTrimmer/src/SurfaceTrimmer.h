#ifdef SURFACETRIMMER_EXPORTS
#define SURFACETRIMMER_API __declspec(dllexport)
#else
#define SURFACETRIMMER_API __declspec(dllimport)
#endif
#pragma once

#include <vector>
#include "MeshData.h"

class SURFACETRIMMER_API SurfaceTrimmer {
public:
	SurfaceTrimmer(void) {};

	// Define the surface trimmer options
	struct Options
	{
		Options() {}
		Options(const float& trim_value = 7, const unsigned int& smooth = 5, const float& aRatio = 0.001) 
			: trim_value(trim_value), smooth(smooth), aRatio(aRatio) {}
		//Parameters
		float trim_value;
		unsigned int smooth;
		float aRatio;
		//To use polygon mesh the AdaptativeSolvers::Mesh structure should change
		//bool polygonMesh;
	};

	bool compute(AdaptativeSolvers::Mesh<float>& mesh_in_out, const Options& options);
};
