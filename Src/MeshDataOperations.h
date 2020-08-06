#pragma once
#include "MeshData.h"
#include <tuple>
#include "PointStreamData.h"

template< typename Real, unsigned int Dim >
struct VertexDataExtractor
{
	typedef PlyVertexWithData< Real, Dim, MultiPointStreamData< Real, PointStreamNormal< Real, Dim >, PointStreamValue< Real >, MultiPointStreamData< Real, PointStreamColor< Real > > > > PlyVertexWithNormalValueAndColor;
	typedef PlyVertexWithData< Real, Dim, MultiPointStreamData< Real, PointStreamNormal< Real, Dim >, PointStreamValue< Real >, MultiPointStreamData< Real > > > PlyVertexWithNormalAndValue;
	typedef PlyVertexWithData< Real, Dim, MultiPointStreamData< Real, PointStreamNormal< Real, Dim >, MultiPointStreamData< Real, PointStreamColor< Real > > > > PlyVertexWithNormalAndColor;
	typedef PlyVertexWithData< Real, Dim, MultiPointStreamData< Real, PointStreamNormal< Real, Dim >, MultiPointStreamData< Real > > > PlyVertexWithNormal;

	template< typename Vertex >
	static void Extract(const Vertex &v, AdaptativeSolvers::Point<Real> &point)
	{
		ERROR_OUT("Unrecognized vertex type");
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
			//Value
			point.value = v.data.template data<1>();
		}
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
			//Value
			point.value = v.data.template data<1>();
		}
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
	typedef PlyVertexWithData< Real, Dim, MultiPointStreamData< Real, PointStreamNormal< Real, Dim >, MultiPointStreamData< Real, PointStreamColor< Real > > > > PlyVertexWithNormalAndColor;
	typedef PlyVertexWithData< Real, Dim, MultiPointStreamData< Real, PointStreamNormal< Real, Dim >, MultiPointStreamData< Real > > > PlyVertexWithNormal;

	template< typename Vertex >
	static void Set(Vertex &v, AdaptativeSolvers::Point<Real> &point)
	{
		ERROR_OUT("Unrecognized vertex type");
	}
	template<>
	static void Set(PlyVertexWithNormalAndColor &v, AdaptativeSolvers::Point<Real> &point)
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
	static void Set(PlyVertexWithNormal &v, AdaptativeSolvers::Point<Real> &point)
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
		for (const auto& point : point_cloud.points)
		{
			Vertex v;
			VertexDataSetter<Real, Dim>::Set(v, *point);
			std::pair< Point< Real, Dim >, TotalPointSampleData > p;
			std::get<0>(p) = v.point;
			std::get<1>(p) = v.data;
			inCorePoints.emplce_back(p);
		}
	}

	template< typename Vertex, typename Real, typename Index>
	static void SetMesh(const AdaptativeSolvers::Mesh<Real>& mesh,
		std::vector< Vertex >& vertices, std::vector< std::vector< Index > >& polygons)
	{
		for (auto point : mesh.points)
		{
			Vertex v;
			VertexDataSetter<Real, 3>::Set(v, point);
			vertices.emplace_back(v);
		}
		for (const auto& face : mesh.faces)
		{
			std::vector< Index > indexes;
			indexes.reserve(3);
			for (size_t i = 0; i < 3; i++)
			{
				indexes.emplace_back(face.point_indices[i]);
			}
			polygons.emplace_back(indexes);
		}
	}

	template< typename Vertex, typename Real, typename Index>
	static void GetMesh(const std::vector< Vertex >& vertices, const std::vector< std::vector< Index > >& polygons,
		AdaptativeSolvers::Mesh<Real>& mesh)
	{
		//clear mesh
		mesh.points.clear();
		mesh.faces.clear();
		mesh.points.reserve(vertices.size());
		mesh.faces.reserve(polygons.size());
		for (auto vert : vertices)
		{
			AdaptativeSolvers::Point<Real> point;
			VertexDataExtractor<Real, 3>::Extract(vert, point);
			mesh.points.emplace_back(point);
		}
		for (auto poly : polygons)
		{
			AdaptativeSolvers::Face face;
			for (size_t i = 0; i < 3; i++)
			{
				face.point_indices[i];
			}
			mesh.faces.emplace_back(face);
		}
	}

	template< typename Vertex, typename Real, typename node_index_type>
	static void GetMesh(CoredMeshData< Vertex, node_index_type >& mesh, AdaptativeSolvers::Mesh<Real>& mesh_in_out)
	{
		mesh->resetIterator();
		typename Vertex::Transform _xForm(iXForm);
		mesh->resetIterator();
		const auto size_in_core_points = mesh->inCorePoints.size();
		const auto size_out_core_points = mesh->outOfCorePointCount();
		//Clear input
		mesh_in_out.points.clear();
		mesh_in_out.faces.clear();
		mesh_in_out.points.reserve(size_in_core_points + size_out_core_points);
		for (size_t i = 0; i < sizeInCorePoints; i++)
		{
			Vertex v = _xForm(mesh->inCorePoints[i]);
			AdaptativeSolvers::Point<Real> point;
			VertexDataExtractor<Real, Dim>::Extract(v, point);
			mesh_in_out.points.emplace_back(point);
		}
		for (size_t i = 0; i < sizeOutOfCorePoints; i++)
		{
			Vertex v;
			mesh->nextOutOfCorePoint(v);
			v = _xForm(v);
			AdaptativeSolvers::Point<Real> point;
			VertexDataExtractor<Real, Dim>::Extract(v, point);
			mesh_in_out.points.emplace_back(point);
		}

		//Write faces
		std::vector< CoredVertexIndex< node_index_type > > polygon;
		const auto size_faces = mesh->polygonCount();
		unsigned int index_vertex;
		mesh_in_out.faces.reserve(size_faces);
		for (size_t i = 0; i < size_faces; i++)
		{
			AdaptativeSolvers::Face face;
			mesh->nextPolygon(polygon);
			const auto size_polygon = polygon.size();
			for (size_t j = 0; j < size_polygon; j++)
			{
				if (polygon[j].inCore)
				{
					index_vertex = polygon[j].idx;
				}
				else
				{
					index_vertex = polygon[j].idx + mesh->inCorePoints.size();
				}
				//Add the pointer to the vertex inside the face
				face.point_indices[j] = index_vertex;
			}
			mesh_in_out.faces.emplace_back(face);
		}
	}
};

