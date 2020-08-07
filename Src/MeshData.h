#pragma once

namespace AdaptativeSolvers
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
		std::vector< unsigned int > point_indices;
	};

	template <typename T>
	struct Mesh
	{
		std::vector< Point<T> > points;
		std::vector< Face > faces;
		//Parameters
		bool has_normal = false;
		bool has_color = false;
		bool has_value = false;
	};
}