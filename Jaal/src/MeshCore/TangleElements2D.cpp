#include "TangleElements2D.hpp"

#include <vector>
#include "Mesh.hpp"
#include "basic_math.hpp"

std::multimap<int,int> boxPairs;
std::map<int,int>  boxFaceID;

int checkPolygon( const vector<Point2D> &polyPnts)
{
    ////////////////////////////////////////////////////////////////////////////
    // This module checks the quality of intersecting polygon between two faces.
    //
    // Since the two input elements are convex ( strict condition): Intersection
    // must be convex and simple. If it is not, then we consider it to be fatal
    // error.
    //
    // If the area of intersection is zero, we return this as an innocuous
    // and empty the contain before returning...
    //
    // The interecting polynomial must be simple. We check convexity first and
    // then simplicity...
    ////////////////////////////////////////////////////////////////////////////

    int npoints = polyPnts.size();
    if( npoints < 3) return -1;

    Polygon_2 poly;
    for( int i = 0; i < npoints; i++) {
        poly.push_back (Point_2 (polyPnts[i][0], polyPnts[i][1]));
    }

    double area = JGeometry::getSignedArea( polyPnts);
    if( fabs(area) < 1.0E-10) return 1;

    if( !poly.is_convex() ) {
        cout << "Warning: Polygon is concave : Area" << area << endl;
    }

    if( !poly.is_simple() ) {
        cout << "Warning: The following polygon is not simple: Results may be incorrect " << endl;
        for( int i = 0; i < npoints; i++) {
            cout << polyPnts[i][0] << " " << polyPnts[i][1] << endl;
        }
        return 1;
    }

    if( poly.orientation() != CGAL::COUNTERCLOCKWISE)
        cout << "Watning: the polygon was not anti-clockwise " << endl;

    return 0;
}

void getPoints( const JFacePtr &face, vector<Point2D> &points)
{
    //
    // Get the points of the face in counter clockwise direction. If the mesh
    // is tangled then some of the elements will have negative orientation,
    // and having positive orientation is prerequisite for overlapping detection.
    //
    short int fsign;
    int err = face->getAttribute("Orient", fsign);
    if( err ) {
        fsign  = JFaceGeometry::getOrientation2D(face);
        face->setAttribute("Orient", fsign);
    }

    int np = face->getSize(0);
    points.resize(np);

    if( fsign > 0) {
        for( int i = 0; i < np; i++) {
            const Point3D &xyz = face->getNodeAt(i)->getXYZCoords();
            points[i][0] = xyz[0];
            points[i][1] = xyz[1];
        }
        return;
    }

    for( int i = 0; i < np; i++) {
        const Point3D &xyz = face->getNodeAt(np-i-1)->getXYZCoords();
        points[i][0] = xyz[0];
        points[i][1] = xyz[1];
    }

    return;
}

/////////////////////////////////////////////////////////////////////////////////////////

bool isAcceptable( const JFacePtr &face)
{
    Polygon_2 poly;

    int npoints = face->getSize(0);
    for( int i = 0; i < npoints; i++) {
        const Point3D xyz = face->getNodeAt(i)->getXYZCoords();
        poly.push_back (Point_2 (xyz[0], xyz[1]));
    }

    // 1St precondition:  The area of the face must not be zero.

    double area = CGAL::to_double( poly.area() );
    if( fabs(area) < 1.0E-10) {
        cout << "Fatal error: Area of a face is closer to zero" << endl;
        exit(0);
    }

    if( !poly.is_simple() ) {
        cout << "Fatal error: The following face is non-simple" << endl;
        for( int i = 0; i < npoints; i++) {
            cout << poly[i].x() << " " << poly[i].y() << endl;
        }
    }

    return 1;
}

////////////////////////////////////////////////////////////////////////////////////////////

int isConvex(const JFacePtr &face)
{
    Polygon_2 poly;

    int npoints = face->getSize(0);
    if( npoints < 4 ) return -1;

    vector<Point3D> qpoints(4);
    for( int i = 0; i < npoints; i++) {
        qpoints[i] = face->getNodeAt(i)->getXYZCoords();
        poly.push_back (Point_2(qpoints[i][0],  qpoints[i][1]));
    }

    if( poly.is_convex() ) return 1;
    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////

int getConcaveCorner(const JFacePtr &face)
{
    int np = face->getSize(0);

    vector<Point3D> qpoints(np);
    for( int i = 0; i < np; i++) {
        qpoints[i] = face->getNodeAt(i)->getXYZCoords();
    }

    Polygon_2 tripoints;
    for( int i = 0; i < np; i++) {
        tripoints.clear();
        tripoints.push_back( Point_2( qpoints[i][0], qpoints[i][1] ));
        tripoints.push_back( Point_2( qpoints[(i+1)%np][0], qpoints[(i+1)%np][1] ));
        tripoints.push_back( Point_2( qpoints[(i+2)%np][0], qpoints[(i+2)%np][1] ));
        double area = CGAL::to_double( tripoints.area() );
        if( area < 0.0) return (i+1)%np;
    }

    return -1;
}

/////////////////////////////////////////////////////////////////////////////////////////
void checkFace( const JFacePtr &face)
{
    Polygon_2 poly;

    int npoints = face->getSize(0);
    for( int i = 0; i < npoints; i++) {
        const Point3D xyz = face->getNodeAt(i)->getXYZCoords();
        poly.push_back (Point_2 (xyz[0], xyz[1]));
    }

    // 1St precondition:  The area of the face must not be zero.

    double area = CGAL::to_double( poly.area() );
    if( fabs(area) < 1.0E-10) {
        cout << "Fatal error: Area of a face is closer to zero" << endl;
        exit(0);
    }

    // We are not supporting self intersecting elements, so every
    // element must be convex...

    if( !poly.is_convex() ) {
        cout << "Fatal Error: Concave elements not supported yet" << endl;
        exit(0);
    }

    // Although, it may be redundant check, because a convex element is
    // simple too, but since I have not considered all the possible cases,
    // I have included this test in the list.

    if( !poly.is_simple() ) {
        cout << "Fatal error: The following face is non-simple" << endl;
        for( int i = 0; i < npoints; i++) {
            cout << poly[i].x() << " " << poly[i].y() << endl;
        }
        exit(0);
    }
}

void BoxBoxIntersectCallback( const Box& a, const Box& b )
{
    int aid = boxFaceID[a.id()];
    int bid = boxFaceID[b.id()];

    int minid = min(aid, bid );
    int maxid = max(aid, bid );

    if( minid < maxid )
        boxPairs.insert(std::pair<int,int>(minid, maxid));
}

/////////////////////////////////////////////////////////////////////////////////////////

template<class Kernel, class Container>
void getPolyPoints(const CGAL::Polygon_with_holes_2<Kernel, Container> &pwh,
                   vector<Point2D> &points)
{
    //
    // Collect the points on the polygon from the datastructure of CGAL...
    //
    typename CGAL::Polygon_2<Kernel, Container> poly;
    typename CGAL::Polygon_2<Kernel, Container>::Vertex_const_iterator vit;

    Point2D xy;
    points.clear();

    if (!pwh.is_unbounded())  {
        poly = pwh.outer_boundary();
        if( poly.size() ) {
            points.reserve(poly.size() );
            for (vit = poly.vertices_begin(); vit != poly.vertices_end(); ++vit) {
                Point_2 pnt = *vit;
                xy[0] = CGAL::to_double(pnt.x());
                xy[1] = CGAL::to_double(pnt.y());
                points.push_back( xy);
            }
        }
    }
    if( points.empty() ) return;

    if( fabs(JGeometry::getSignedArea( points)) <= 1.0E-12) {
        points.clear();
        return;
    }

    // This is just a sanity check, because it will always pass through, because
    // the intersection of two convex shape is always a convex region..

    Polygon_2 P;
    int nsize = points.size();
    for( int i = 0; i < nsize; i++)
        P.push_back (Point_2 (points[i][0], points[i][1]));

    if( P.is_convex() ) {
        vector<Point2D>::iterator it;
        it = std::unique( points.begin(), points.end() );
        points.erase(it, points.end() );
        if( points.back() == points.front() )
            points.pop_back();
    }

}

/////////////////////////////////////////////////////////////////////////////////////////

bool hasIntersection( const vector<Point2D> &aPoints, const vector<Point2D> &bPoints)
{
    //
    // Check whether tiangle A and triangle B overlap, if they do, then return the
    // points of the polygon formed from the intersection ...
    // The points must be counter-clockwisely oriented. The outpoint points are in the
    // counter clockwise direction ....
    //
    int nsize;

    Polygon_2 P;
    nsize = aPoints.size();
    for( int i = 0; i < nsize; i++)
        P.push_back (Point_2 (aPoints[i][0], aPoints[i][1]));

    Polygon_2 Q;
    nsize = bPoints.size();
    for( int i = 0; i < nsize; i++)
        Q.push_back (Point_2 (bPoints[i][0], bPoints[i][1]));

    Polylist polylist;
    return CGAL::do_intersect(P, Q);
}

/////////////////////////////////////////////////////////////////////////////////////////

void getIntersection( const vector<Point2D> &aPoints, const vector<Point2D> &bPoints,
                      vector<Point2D> &cPoints)
{
    // Check whether tiangle A and triangle B overlap, if they do, then return the
    // points of the polygon formed from the intersection ...
    // The points must be counter-clockwisely oriented. The outpoint points are in the
    // counter clockwise direction ....
    //
    cPoints.clear();

    int nsize;

    Polygon_2 P;
    nsize = aPoints.size();
    for( int i = 0; i < nsize; i++)
        P.push_back (Point_2 (aPoints[i][0], aPoints[i][1]));

    Polygon_2 Q;
    nsize = bPoints.size();
    for( int i = 0; i < nsize; i++)
        Q.push_back (Point_2 (bPoints[i][0], bPoints[i][1]));

    if( CGAL::do_intersect(P, Q) ) {
        Polylist polylist;
        CGAL::intersection (Q, P, std::back_inserter(polylist));

        Polylist::const_iterator  it;
        for (it = polylist.begin(); it != polylist.end(); ++it)
            getPolyPoints(*it, cPoints);

        checkPolygon( cPoints);
    }
}

/////////////////////////////////////////////////////////////////////////////////////////

void getIntersection( const JFacePtr &face1, const JFacePtr &face2, vector<Point2D> &polyPoints)
{
    polyPoints.clear();
    assert( face1 != face2);

    vector<Point2D> aPoints, bPoints;
    getPoints( face1, aPoints);
    getPoints( face2, bPoints);

    getIntersection( aPoints, bPoints, polyPoints);
}
/////////////////////////////////////////////////////////////////////////////////////////

int isInside( const vector<Point2D> &polyPoints, const Point2D &qPoint)
{
    vector<Point_2> points;
    int np = polyPoints.size();

    points.resize(np);
    for( int i = 0; i < np; i++)
        points[i] = Point_2(polyPoints[i][0], polyPoints[i][1]);

    Point_2 testPoint( qPoint[0], qPoint[1] );

    int stat;
    stat = CGAL::bounded_side_2( points.begin(), points.end(), testPoint, Kernel() );

    if( stat == CGAL::ON_UNBOUNDED_SIDE) {
        cout << setprecision(10) << fixed;
        cout << "Point is outside the polygon" << endl;
        for( int i = 0; i < np; i++)
            cout << polyPoints[i][0] << "  " << polyPoints[i][1] << endl;
        cout << "Query Point " << endl;
        cout << testPoint[0] << "  " << testPoint[1] << endl;
        exit(0);
        return 0;
    }

    return 1;
}
//
/////////////////////////////////////////////////////////////////////////////////////////

int isInside( const JFacePtr &face, const Point2D &qPoint)
{
    vector<Point2D> points;
    getPoints( face, points);
    return isInside( points, qPoint);
}
/////////////////////////////////////////////////////////////////////////////////////////


