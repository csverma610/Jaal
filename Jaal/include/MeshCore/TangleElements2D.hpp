#pragma once

#include <vector>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/Bbox_2.h>

#include "Mesh.hpp"

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef Kernel::Point_2                     Point_2;
typedef CGAL::Polygon_2<Kernel>             Polygon_2;
typedef CGAL::Polygon_with_holes_2<Kernel>  Polygon_with_holes_2;
typedef std::list<Polygon_with_holes_2>     Polylist;

typedef CGAL::Box_intersection_d::Box_d<double,2> Box;
typedef CGAL::Bbox_2  Bbox;

int  isConvex(const JFacePtr &face);
int  isInside( const JFacePtr &face, const Point2D &qPoint);
int  isInside( const vector<Point2D> &polyPoints, const Point2D &qPoint);
void getIntersection( const JFacePtr &face1, const JFacePtr &face2, vector<Point2D> &polyPoints);
void checkFace( const JFacePtr &f);

/*
#include "MeshEntity.hpp"
using namespace std;

int  checkPolygon( const vector<Point2D> &polyPnts);
int  getConcaveCorner(const JFacePtr &face);
bool isAcceptable( const JFacePtr &face);
void BoxBoxIntersectCallback( const Box& a, const Box& b);
bool hasIntersection( const vector<Point2D> &aPoints, const vector<Point2D> &bPoints);
void getIntersection( const vector<Point2D> &aPoints, const vector<Point2D> &bPoints,
                      vector<Point2D> &cPoints);
*/
/////////////////////////////////////////////////////////////////////////////////////////
