#pragma once
#include "MeshData.h"
#include <tuple>
#include "PointStreamData.h"

template< typename Real, unsigned int Dim >
struct VertexDataExtractor
{
	typedef PlyVertexWithData< Real, Dim, MultiPointStreamData< Real, PointStreamValue< Real >, PointStreamNormal< Real, Dim >, PointStreamColor< float > > > PlyVertexWithValueNormalAndColor;
	typedef PlyVertexWithData< Real, Dim, MultiPointStreamData< Real, PointStreamValue< Real >, PointStreamNormal< Real, Dim > > > PlyVertexWithValueAndNormal;
	typedef PlyVertexWithData< Real, Dim, MultiPointStreamData< Real, PointStreamValue< Real >, PointStreamColor< float > > > PlyVertexWithValueAndColor;
	typedef PlyVertexWithData< Real, Dim, MultiPointStreamData< Real, PointStreamValue< Real > > > PlyVertexWithValue;

	typedef PlyVertexWithData< Real, Dim, MultiPointStreamData< Real, PointStreamNormal< Real, Dim >, PointStreamValue< Real >, MultiPointStreamData< Real, PointStreamColor< float > > > > PlyVertexWithNormalValueAndColor;
	typedef PlyVertexWithData< Real, Dim, MultiPointStreamData< Real, PointStreamNormal< Real, Dim >, PointStreamValue< Real >, MultiPointStreamData< Real > > > PlyVertexWithNormalAndValue;
	typedef PlyVertexWithData< Real, Dim, MultiPointStreamData< Real, PointStreamNormal< Real, Dim >, MultiPointStreamData< Real, PointStreamColor< Real > > > > PlyVertexWithNormalAndColor;
	typedef PlyVertexWithData< Real, Dim, MultiPointStreamData< Real, PointStreamNormal< Real, Dim >, MultiPointStreamData< Real > > > PlyVertexWithNormal;

	template< typename Vertex >
	static void Extract(const Vertex &v, AdaptativeSolvers::Point<Real> &point)
	{
		ERROR_OUT("Unrecognized vertex type");
	}
	template<>
	static void Extract(const PlyVertexWithValueNormalAndColor &v, AdaptativeSolvers::Point<Real> &point)
	{
		for (size_t i = 0; i < 3; i++)
		{
			//Point
			point.xyz[i] = v.point.coords[i];
			//Normal
			point.normal[i] = v.data.template data<1>().coords[i];
			//Color 
			point.color[i] = v.data.template data<2>().coords[i];
		}
		//Value
		point.value = v.data.template data<0>();
	}
	template<>
	static void Extract(const PlyVertexWithValueAndNormal &v, AdaptativeSolvers::Point<Real> &point)
	{
		for (size_t i = 0; i < 3; i++)
		{
			//Point
			point.xyz[i] = v.point.coords[i];
			//Normal
			point.normal[i] = v.data.template data<1>().coords[i];
		}
		//Value
		point.value = v.data.template data<0>();
	}
	template<>
	static void Extract(const PlyVertexWithValueAndColor &v, AdaptativeSolvers::Point<Real> &point)
	{
		for (size_t i = 0; i < 3; i++)
		{
			//Point
			point.xyz[i] = v.point.coords[i];
			//Color
			point.color[i] = v.data.template data<1>().coords[i];
		}
		//Value
		point.value = v.data.template data<0>();
	}
	template<>
	static void Extract(const PlyVertexWithValue &v, AdaptativeSolvers::Point<Real> &point)
	{
		for (size_t i = 0; i < 3; i++)
		{
			//Point
			point.xyz[i] = v.point.coords[i];
		}
		//Value
		point.value = v.data.template data<0>();
	}
	template<>
	static void Extract(const PlyVertexWithNormalValueAndColor &v, AdaptativeSolvers::Point<Real> &point)
	{
		for (size_t i = 0; i < 3; i++)
		{
			//Point
			point.xyz[i] = v.point.coords[i];
			//Normal
			point.normal[i] = v.data.template data<0>().coords[i];
			//Color 
			point.color[i] = std::get< 0 >(v.data.template data<2>()).data().coords[i];
		}
		//Value
		point.value = v.data.template data<1>();
	}
	template<>
	static void Extract(const PlyVertexWithNormalAndValue &v, AdaptativeSolvers::Point<Real> &point)
	{
		for (size_t i = 0; i < 3; i++)
		{
			//Point
			point.xyz[i] = v.point.coords[i];
			//Normal
			point.normal[i] = v.data.template data<0>().coords[i];
		}
		//Value
		point.value = v.data.template data<1>();
	}
	template<>
	static void Extract(const PlyVertexWithNormalAndColor &v, AdaptativeSolvers::Point<Real> &point)
	{
		for (size_t i = 0; i < 3; i++)
		{
			//Point
			point.xyz[i] = v.point.coords[i];
			//Normal
			point.normal[i] = v.data.template data<0>().coords[i];
			//Color 
			point.color[i] = std::get< 0 >(v.data.template data<1>()).data().coords[i];
		}
	}
	template<>
	static void Extract(const PlyVertexWithNormal &v, AdaptativeSolvers::Point<Real> &point)
	{
		for (size_t i = 0; i < 3; i++)
		{
			//Point
			point.xyz[i] = v.point.coords[i];
			//Normal
			point.normal[i] = v.data.template data<0>().coords[i];
		}
	}
};

template< typename Real, unsigned int Dim >
struct VertexDataSetter
{
	typedef PlyVertexWithData< Real, Dim, MultiPointStreamData< Real, PointStreamValue< Real >, PointStreamNormal< Real, Dim >, PointStreamColor< float > > > PlyVertexWithValueNormalAndColor;
	typedef PlyVertexWithData< Real, Dim, MultiPointStreamData< Real, PointStreamValue< Real >, PointStreamNormal< Real, Dim > > > PlyVertexWithValueAndNormal;
	typedef PlyVertexWithData< Real, Dim, MultiPointStreamData< Real, PointStreamValue< Real >, PointStreamColor< float > > > PlyVertexWithValueAndColor;
	typedef PlyVertexWithData< Real, Dim, MultiPointStreamData< Real, PointStreamValue< Real > > > PlyVertexWithValue;

	typedef PlyVertexWithData< Real, Dim, MultiPointStreamData< Real, PointStreamNormal< Real, Dim >, PointStreamValue< Real >, MultiPointStreamData< Real, PointStreamColor< float > > > > PlyVertexWithNormalValueAndColor;
	typedef PlyVertexWithData< Real, Dim, MultiPointStreamData< Real, PointStreamNormal< Real, Dim >, PointStreamValue< Real >, MultiPointStreamData< Real > > > PlyVertexWithNormalAndValue;
	typedef PlyVertexWithData< Real, Dim, MultiPointStreamData< Real, PointStreamNormal< Real, Dim >, MultiPointStreamData< Real, PointStreamColor< float > > > > PlyVertexWithNormalAndColor;
	typedef PlyVertexWithData< Real, Dim, MultiPointStreamData< Real, PointStreamNormal< Real, Dim >, MultiPointStreamData< Real > > > PlyVertexWithNormal;

	template< typename Vertex >
	static void Set(Vertex &v, const AdaptativeSolvers::Point<Real> &point)
	{
		ERROR_OUT("Unrecognized vertex type");
	}
	template<>
	static void Set(PlyVertexWithValueNormalAndColor &v, const AdaptativeSolvers::Point<Real> &point)
	{
		for (size_t i = 0; i < 3; i++)
		{
			//Point
			v.point.coords[i] = point.xyz[i];
			//Normal
			v.data.template data<1>().coords[i] = point.normal[i];
			//Color 
			v.data.template data<2>().coords[i] = point.color[i];
		}
		//Value
		v.data.template data<0>() = point.value;
	}
	template<>
	static void Set(PlyVertexWithValueAndNormal &v, const AdaptativeSolvers::Point<Real> &point)
	{
		for (size_t i = 0; i < 3; i++)
		{
			//Point
			v.point.coords[i] = point.xyz[i];
			//Normal
			v.data.template data<1>().coords[i] = point.normal[i];
		}
		//Value
		v.data.template data<0>() = point.value;
	}
	template<>
	static void Set(PlyVertexWithValueAndColor &v, const AdaptativeSolvers::Point<Real> &point)
	{
		for (size_t i = 0; i < 3; i++)
		{
			//Point
			v.point.coords[i] = point.xyz[i];
			//Color
			v.data.template data<1>().coords[i] = point.color[i];
		}
		//Value
		v.data.template data<0>() = point.value;
	}
	template<>
	static void Set(PlyVertexWithValue &v, const AdaptativeSolvers::Point<Real> &point)
	{
		for (size_t i = 0; i < 3; i++)
		{
			//Point
			v.point.coords[i] = point.xyz[i];
		}
		//Value
		v.data.template data<0>() = point.value;
	}
	template<>
	static void Set(PlyVertexWithNormalValueAndColor &v, const AdaptativeSolvers::Point<Real> &point)
	{
		for (size_t i = 0; i < 3; i++)
		{
			//Point
			v.point.coords[i] = point.xyz[i];
			//Normal
			v.data.template data<0>().coords[i] = point.normal[i];
			//Color 
			std::get< 0 >(v.data.template data<2>()).data().coords[i] = point.color[i];
		}
		//Value
		v.data.template data<1>() = point.value;
	}
	template<>
	static void Set(PlyVertexWithNormalAndValue &v, const AdaptativeSolvers::Point<Real> &point)
	{
		for (size_t i = 0; i < 3; i++)
		{
			//Point
			v.point.coords[i] = point.xyz[i];
			//Normal
			v.data.template data<0>().coords[i] = point.normal[i];
		}
		//Value
		v.data.template data<1>() = point.value;
	}
	template<>
	static void Set(PlyVertexWithNormalAndColor &v, const AdaptativeSolvers::Point<Real> &point)
	{
		for (size_t i = 0; i < 3; i++)
		{
			//Point
			v.point.coords[i] = point.xyz[i];
			//Normal
			v.data.template data<0>().coords[i] = point.normal[i];
			//Color 
			std::get< 0 >(v.data.template data<1>()).data().coords[i] = point.color[i];
		}
	}
	template<>
	static void Set(PlyVertexWithNormal &v, const AdaptativeSolvers::Point<Real> &point)
	{
		for (size_t i = 0; i < 3; i++)
		{
			//Point
			v.point.coords[i] = point.xyz[i];
			//Normal
			v.data.template data<0>().coords[i] = point.normal[i];
		}
	}
};


class SetGetMesh
{
public:
	SetGetMesh() {};
	~SetGetMesh() {};

	template< typename Vertex, typename Real, unsigned int Dim, typename TotalPointSampleData>
	static void SetPoints(const AdaptativeSolvers::Mesh<Real>& point_cloud, std::vector< std::pair< Point< Real, Dim >, TotalPointSampleData > >& inCorePoints)
	{
		const auto points_size = point_cloud.points.size();
		inCorePoints.resize(points_size);
#pragma omp parallel for
		for (int i = 0; i < points_size; i++)
		{
			Vertex v;
			VertexDataSetter<Real, Dim>::Set(v, point_cloud.points[i]);
			std::pair< Point< Real, Dim >, TotalPointSampleData > p;
			std::get<0>(p) = v.point;
			std::get<1>(p) = v.data;
			inCorePoints[i] = p;
		}
	}

	template< typename Vertex, typename Real, typename Index>
	static void SetMesh(const AdaptativeSolvers::Mesh<Real>& mesh,
		std::vector< Vertex >& vertices, std::vector< std::vector< Index > >& polygons)
	{
		const auto points_size = mesh.points.size();
		vertices.resize(points_size);
#pragma omp parallel for
		for (int i = 0; i < points_size; i++)
		{
			Vertex v;
			VertexDataSetter<Real, 3>::Set(v, mesh.points[i]);
			vertices[i] = v;
		}
		const auto faces_size = mesh.faces.size();
		polygons.resize(faces_size);
#pragma omp parallel for
		for (int i = 0; i < faces_size; i++)
		{
			std::vector< Index > indexes;
			indexes.reserve(mesh.faces[i].point_indices.size());
			for (auto idx : mesh.faces[i].point_indices)
			{
				indexes.emplace_back(idx);
			}
			polygons[i] = indexes;
		}
	}

	template< typename Vertex, typename Real, typename Index>
	static void GetMesh(const std::vector< Vertex >& vertices, const std::vector< std::vector< Index > >& polygons,
		AdaptativeSolvers::Mesh<Real>& mesh)
	{
		//clear mesh
		mesh.points.clear();
		mesh.faces.clear();
		const auto vertices_size = vertices.size();
		const auto polygons_size = polygons.size();
		mesh.points.resize(vertices_size);
		mesh.faces.resize(polygons_size);
#pragma omp parallel for
		for (int i = 0; i < vertices_size; i++)
		{
			AdaptativeSolvers::Point<Real> point;
			VertexDataExtractor<Real, 3>::Extract(vertices[i], point);
			mesh.points[i] = point;
		}
#pragma omp parallel for
		for (int i = 0; i < polygons_size; i++)
		{
			AdaptativeSolvers::Face face;
			face.point_indices.reserve(polygons[i].size());
			for (auto idx : polygons[i])
			{
				face.point_indices.emplace_back(idx);
			}
			mesh.faces[i] = face;
		}
	}

	template< typename Vertex, typename Real, typename node_index_type, unsigned int Dim, unsigned int ... FEMSigs>
	static void GetMesh(CoredMeshData< Vertex, node_index_type >& mesh, AdaptativeSolvers::Mesh<Real>& mesh_in_out, XForm< Real, sizeof...(FEMSigs) + 1 > unitCubeToModel)
	{
		typename Vertex::Transform _xForm(unitCubeToModel);
		mesh.resetIterator();
		const auto size_in_core_points = mesh.inCorePoints.size();
		const auto size_out_core_points = mesh.outOfCorePointCount();
		//Clear input
		mesh_in_out.points.clear();
		mesh_in_out.faces.clear();
		mesh_in_out.points.resize(size_in_core_points + size_out_core_points);
#pragma omp parallel for
		for (int i = 0; i < size_in_core_points; i++)
		{
			Vertex v = _xForm(mesh.inCorePoints[i]);
			AdaptativeSolvers::Point<Real> point;
			VertexDataExtractor<Real, Dim>::Extract(v, point);
			mesh_in_out.points[i] = point;
		}
		for (size_t i = 0; i < size_out_core_points; i++)
		{
			Vertex v;
			mesh.nextOutOfCorePoint(v);
			v = _xForm(v);
			AdaptativeSolvers::Point<Real> point;
			VertexDataExtractor<Real, Dim>::Extract(v, point);
			mesh_in_out.points[i] = point;
		}

		//Write faces
		std::vector< CoredVertexIndex< node_index_type > > polygon;
		const auto size_faces = mesh.polygonCount();
		unsigned int index_vertex;
		mesh_in_out.faces.reserve(size_faces);
		for (size_t i = 0; i < size_faces; i++)
		{
			AdaptativeSolvers::Face face;
			mesh.nextPolygon(polygon);
			const auto size_polygon = polygon.size();
			face.point_indices.reserve(size_polygon);
			for (size_t j = 0; j < size_polygon; j++)
			{
				if (polygon[j].inCore)
				{
					index_vertex = polygon[j].idx;
				}
				else
				{
					index_vertex = polygon[j].idx + mesh.inCorePoints.size();
				}
				//Add the pointer to the vertex inside the face
				face.point_indices.emplace_back(index_vertex);
			}
			mesh_in_out.faces.emplace_back(face);
		}
	}
};

