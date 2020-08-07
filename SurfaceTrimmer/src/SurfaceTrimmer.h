#ifdef SURFACETRIMMER_EXPORTS
#define SURFACETRIMMER_API __declspec(dllexport)
#else
#define SURFACETRIMMER_API __declspec(dllimport)
#endif
#pragma once

#include <vector>
#include "MeshData.h"

class SURFACETRIMMER_API SurfaceTrimmer 
{
public:
	SurfaceTrimmer(void) {};

	// Define the surface trimmer options
	struct Options
	{
		Options(const float& trim_value = 7, const unsigned int& smooth = 5, const float& aRatio = 0.001, const bool& polygonMesh = false)
			: trim_value(trim_value), smooth(smooth), aRatio(aRatio), polygonMesh(polygonMesh) {}
		//Parameters
		float trim_value;
		unsigned int smooth;
		float aRatio;
		bool polygonMesh;
	};

	bool compute(AdaptativeSolvers::Mesh<float>& mesh_in_out, const Options& options);
};
