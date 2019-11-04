#include <iomanip>
#include "FEMSpace.hpp"

///////////////////////////////////////////////////////////////////////////////
int JEdgeElement :: setGaussPoints(int n)
{
    int ngauss;
    switch( order ) {
    case 1:
        ngauss = max(1,n);
        break;
    case 2:
        ngauss = max(3,n);
        break;
    }

    gaussPoints.clear();

    vector<double> u, w;
    u.reserve(ngauss);
    w.reserve(ngauss);

    switch( ngauss ) {
    case 1:
        u += 0.0;
        w += 2.0;
        break;
    case 2:
        u += -0.577350269189626, 0.577350269189626;
        w += 1.0, 1.0;
        break;
    case 3:
        u += -0.774596669241483, 0.0,  0.774596669241483;
        w +=  0.555555555555556, 0.88888888888888889, 0.5555555555555556;
        break;
    case 4:
        u += -0.8611363115, -0.3399810435, 0.3399810435, 0.8611363115;
        w += 0.3478548451, 0.6521451548, 0.6521451548, 0.3478548451;
        break;
    case 5:
        u += -0.906179846, -0.53846931, 0,  0.53846931, 0.906179846;
        w += 0.236926885, 0.4786286704, 0.5688888888, 0.4786286704, 0.236926885;
        break;
    default:
        cout << "Fatal Error: Invalud number of Gauss points on ed=ge " << endl;
        exit(0);
    }

    gaussPoints.resize(ngauss);
    for( int i = 0; i < ngauss; i++) {
        gaussPoints[i][0] = u[i];
        gaussPoints[i][1] = w[i];
    }

#ifdef DEBUG
    double sum = 0.0;
    for( int i = 0; i < ngauss; i++)
        sum += w[i];
    assert( fabs(sum - 2.0) < 1.0E-05) );
#endif
    return 0;
}

////////////////////////////////////////////////////////////////////////////

int JEdgeElement :: getShapeFunc(double u, vector<double> &N)
{
    N.clear();
    switch( order) {
    case 1:
        N.resize(2);
        N[0]  = 0.5*(1.0-u);
        N[1]  = 0.5*(1.0+u);
        return 0;
    case 2:
        N.resize(3);
        N[0] = 0.5*u*(u-1.0);
        N[1] = (1.0-u)*(1.0+u);
        N[2] = 0.5*u*(u+1.0);
        return 0;
    }
    cout << "Fatal Error: Shape function not evaluate " << endl;
    return 1;
}

////////////////////////////////////////////////////////////////////////////

int JEdgeElement:: getShapeDeriv(double u, vector<double> &gradN)
{
    gradN.clear();

    switch( order ) {
    case 1:
        gradN.resize(2);
        gradN[0] =  -1/2;
        gradN[1] =   1/2;
        return 0;
    case 2:
        gradN.resize(3);
        gradN[0] =  0.5*(2*u-1);
        gradN[1] =  (-2*u);
        gradN[2] =  0.5*(2*u+1);
        return 0;
    }

    cout << "Fatal Error: Shape gradient  not evaluate " << endl;
    return 1;
}

////////////////////////////////////////////////////////////////////////////

int JTriElement:: setGaussPoints(int n)
{
    int ngauss;
    switch( order ) {
    case 1:
        ngauss  = max(1,n);
        break;
    case 2:
        ngauss  = max(3,n);
        break;
    }
    cout << "Order " << order << endl;

    vector<double> u, v, w;

    u.reserve(ngauss);
    v.reserve(ngauss);
    w.reserve(ngauss);

    switch( ngauss ) {
    case 1:
        u += 1.0/3.0;
        v += 1.0/3.0;
        w += 1.0/2.0;
        break;
    case 3:
        u += 1.0/6.0, 2.0/3.0, 1.0/6.0;
        v += 1.0/6.0, 1.0/6.0, 2.0/3.0;
        w += 1.0/6.0, 1.0/6.0, 1.0/6.0;
        break;
    case 6:
        u += 0.44594849091597, 0.44594849091597, 0.10810301816807, 0.09157621350977, 0.09157621350977, 0.81684757298046;
        v += 0.44594849091597, 0.10810301816807, 0.44594849091597, 0.09157621350977, 0.81684757298046, 0.09157621350977;
        w += 0.11169079483901, 0.11169079483901, 0.11169079483901, 0.05497587182766, 0.05497587182766, 0.05497587182766;
        break;
    case 7:
        u += 0.1012865073235, 0.7974269853531, 0.1012865073235, 0.4701420641051,
             0.4701420641051, 0.0597158717898, 0.3333333333333;
        v += 0.1012865073235, 0.1012865073235, 0.7974269853531, 0.0597158717898,
             0.4701420641051, 0.4701420641051, 0.3333333333333;
        w += 0.0629695902724, 0.0629695902724, 0.0629695902724, 0.0661970763942,
             0.0661970763942, 0.0661970763942, 0.1125000000000;
        break;
    case 12:
        u += 0.24928674517091, 0.24928674517091, 0.50142650965818, 0.06308901449150,
             0.06308901449150, 0.87382197101700, 0.31035245103378, 0.63650249912140,
             0.05314504984482, 0.63650249912140, 0.31035245103378, 0.05314504984482;
        v += 0.24928674517091, 0.50142650965818, 0.24928674517091, 0.06308901449150,
             0.87382197101700, 0.06308901449150, 0.63650249912140, 0.05314504984482,
             0.31035245103378, 0.31035245103378, 0.05314504984482, 0.63650249912140;
        w += 0.05839313786319, 0.05839313786319, 0.05839313786319, 0.02542245318510,
             0.02542245318510, 0.02542245318510, 0.04142553780919, 0.04142553780919,
             0.04142553780919, 0.04142553780919, 0.04142553780919, 0.04142553780919;
        break;
    default:
        cout << "Warning: Invalid Guass points for triangle element(1,3,7, 13)" << endl;
        return 1;
    }

    gaussPoints.resize(ngauss);
    for( int i = 0; i < ngauss; i++) {
        gaussPoints[i][0] = u[i];
        gaussPoints[i][1] = v[i];
        gaussPoints[i][2] = w[i];
    }

#ifdef DEBUG
    double sum = 0.0;
    for( int i = 0; i < ngauss; i++)
        sum += w[i];
    assert( fabs(sum - 0.5) < 1.0E-05) );
#endif

    return 0;
}

//////////////////////////////////////////////////////////////////////////////////

int JTriElement:: getShapeFunc( const Point2D &uv, vector<double> &N)
{
    double u = uv[0];
    double v = uv[1];

    if( order == 1 ) {
        N.resize(3);
        N[0] = 1- u - v;
        N[1] = u;
        N[2] = v;
        return 0;
    }

    N.resize(6);
    double lambda = 1-u-v;
    N[0] = lambda*(2*lambda-1);
    N[1] = u*(2*u-1);
    N[2] = v*(2*v-1);
    N[3] =  4*u*lambda;
    N[4] =  4*u*v;
    N[5] =  4*v*lambda;
    return 0;
}

//////////////////////////////////////////////////////////////////////////////////

int JTriElement :: getShapeDeriv( const Point2D &uv,  MatrixXd &gradN)
{
    double u  = uv[0];
    double v  = uv[1];
    switch( order ) {
    case 1:
        gradN.resize(2,3);
        gradN(0,0) = -1;
        gradN(0,1) = 1;
        gradN(0,2) = 0;
        gradN(1,0) = -1;
        gradN(1,1) = 0;
        gradN(1,2) = 1;
        return 0;
    case 2:
        gradN.resize(2,6);
        gradN(0,0)  =  -3+4*u+4*v;
        gradN(0,1)  =   4*u-1;
        gradN(0,2)  =   0 ;
        gradN(0,3)  =   4-8*u-4*v ;
        gradN(0,4)  =   4*v;
        gradN(0,5)  =  -4*v;

        gradN(1,0)  =  -3+4*u+4*v;
        gradN(1,1)  =   0.0;
        gradN(1,2)  =  4*v-1;
        gradN(1,3)  = -4*u;
        gradN(1,4)  =  4*u;
        gradN(1,5)  = 4-4*u-8*v;
        return 0;
    }

    cout << "Fatal Error: Shape derivatives not evaluated " << endl;
    return 1;
}
//////////////////////////////////////////////////////////////////////////////////
double JTriElement :: eval( const Point2D &uv, const vector<double> &phi)
{
    vector<double> N;
    int err = getShapeFunc(uv, N);
    if( err ) return 0.0;

    assert( N.size() == phi.size() );

    int numnodes = phi.size();
    double sum = 0.0;
    for( int i = 0; i < numnodes; i++)
        sum += N[i]*phi[i];
    return sum;
}

/////////////////////////////////////////////////////////////////////////////////
Array2D JTriElement :: eval( const Point2D &uv, const vector<Array2D> &phi)
{
    Array2D val;
    vector<double> N;
    int err = getShapeFunc(uv, N);
    if( err ) return val;

    assert( N.size() == phi.size() );

    int numnodes = phi.size();
    double xsum = 0.0;
    double ysum = 0.0;
    for( int i = 0; i < numnodes; i++) {
        xsum += N[i]*phi[i][0];
        xsum += N[i]*phi[i][1];
    }
    val[0] = xsum;
    val[1] = ysum;
    return val;
}
/////////////////////////////////////////////////////////////////////////////////

Array3D JTriElement :: eval( const Point2D &uv, const vector<Array3D> &phi)
{
    Array3D val;
    vector<double> N;
    int err = getShapeFunc(uv, N);
    if( err ) return val;

    assert( N.size() == phi.size() );

    int numnodes = phi.size();
    double xsum = 0.0;
    double ysum = 0.0;
    double zsum = 0.0;
    for( int i = 0; i < numnodes; i++) {
        xsum += N[i]*phi[i][0];
        ysum += N[i]*phi[i][1];
        zsum += N[i]*phi[i][2];
    }
    val[0] = xsum;
    val[1] = ysum;
    val[2] = zsum;
    return val;

}
/////////////////////////////////////////////////////////////////////////////////

int JTriElement:: getXYCoords( const vector<Point2D> &xyCoords, const Point2D &uv, Point2D &xy)
{
    double u = uv[0];
    double v = uv[1];
    double eps = 1.0E-06;

    if( u < 0.0 ||  u > 1.0) return 1;
    if( v < 0.0 ||  v > 1.0) return 2;

    if( xyCoords.size() != 3 ) {
        cout << "Warning: Tri element must have three nodes " << endl;
        return 3;
    }

    vector<double>  N;
    getShapeFunc( uv, N);

    double xsum = 0.0;
    double ysum = 0.0;

    for( int i = 0; i < N.size(); i++) {
        xsum += xyCoords[i][0]*N[i];
        ysum += xyCoords[i][1]*N[i];
    }

    xy[0] = xsum;
    xy[1] = ysum;

    return 0;
}

/////////////////////////////////////////////////////////////////////////////////

int JQuadElement :: setGaussPoints(int n)
{
    int ngauss;
    switch( order ) {
    case 1:
        ngauss = max(4,n);
        break;
    case 2:
        ngauss = max(9,n);
        break;
    }

    vector<double> u, v, w;

    u.reserve(ngauss);
    v.reserve(ngauss);
    w.reserve(ngauss);

    double k;
    switch( ngauss) {
    case 4:
        k = 1/sqrt(3.0);
        u += -k, k, k,-k;
        v += -k,-k, k, k;;
        w +=  1.0, 1.0, 1.0, 1.0;
        break;
    case 9:
        k  = sqrt(3.0/5.0);
        u += -k, 0.0, k, -k, 0.0, k, -k, 0.0, k;
        v += -k, -k, -k, 0.0, 0.0, 0.0, k, k, k;
        w += 25.0/81.0, 40.0/81.0, 25.0/81.0,
             40.0/81.0, 64.0/81.0, 40.0/81.0,
             25.0/81.0, 40.0/81.0, 25.0/81.0;
        break;
    default:
        cout << "Warning: Invalid number of Gauss points for Quad element(4,9)" << endl;
        return 1;
    }

    gaussPoints.resize(ngauss);
    for( int i = 0; i < ngauss; i++) {
        gaussPoints[i][0] = u[i];
        gaussPoints[i][1] = v[i];
        gaussPoints[i][2] = w[i];
    }

#ifdef DEBUG
    double sum = 0.0;
    for( int i = 0; i < ngauss; i++)
        sum += w[i];
    assert( fabs(sum - 4.0) < 1.0E-05);
#endif
    return 0;
}

////////////////////////////////////////////////////////////////////////////////

int JQuadElement :: getShapeFunc( const Point2D &uv, vector<double> &N)
{
    N.clear();
    double u = uv[0];
    double v = uv[1];

    if( order == 1) {
        N.resize(4);
        N[0] = 0.25*(1.0-u)*(1.0-v);
        N[1] = 0.25*(1.0+u)*(1.0-v);
        N[2] = 0.25*(1.0+u)*(1.0+v);
        N[3] = 0.25*(1.0-u)*(1.0+v);
        return 0;
    }

    N.resize(8);
    N[0] = 0.25*(1.0-u)*(v-1.0)*(u+v+1.0);
    N[1] = 0.25*(1.0+u)*(v-1.0)*(v-u+1.0);
    N[2] = 0.25*(1.0+u)*(1.0+v)*(u+v-1.0);
    N[3] = 0.25*(u-1.0)*(v+1.0)*(u-v+1.0);
    N[4] = 0.50*(1.0-v)*(1.0-u*u);
    N[5] = 0.50*(1.0+u)*(1.0-v*v);
    N[6] = 0.50*(1.0+v)*(1.0-u*u);
    N[7] = 0.50*(1.0-u)*(1.0-v*v);
    return 0;
}
////////////////////////////////////////////////////////////////////////////////

int JQuadElement :: getShapeDeriv(const Point2D &uv, MatrixXd  &gradN)
{
    double u = uv[0];
    double v = uv[1];

    if( order == 1 ) {
        gradN.resize(2,4);
        gradN(0,0) = -0.25*(1.0-v);
        gradN(0,1) =  0.25*(1.0-v);
        gradN(0,2) =  0.25*(1.0+v);
        gradN(0,3) = -0.25*(1.0+v);

        gradN(1,0) = -0.25*(1.0-u);
        gradN(1,1) = -0.25*(1.0+u);
        gradN(1,2) =  0.25*(1.0+u);
        gradN(1,3) =  0.25*(1.0-u);
        return 0;
    }

    gradN.resize(2,8);
    gradN(0,0)  = 0.25*(-(u-1.0)*(v-1.0) - (v-1.0)*(u + v + 1.0));
    gradN(1,0)  = 0.25*(-(u-1.0)*(v-1.0) - (u-1.0)*(u + v + 1.0));

    gradN(0,1)  = 0.25*(-(u+1.0)*(v-1.0) - (v-1.0)*(u - v - 1.0));
    gradN(1,1)  = 0.25*(+(u+1.0)*(v-1.0) - (u+1.0)*(u - v - 1.0));

    gradN(0,2)  = 0.25*((u+1.0)*(v+1.0) + (v+1.0)*(u + v - 1.0));
    gradN(1,2)  = 0.25*((u+1.0)*(v+1.0) + (u+1.0)*(u + v - 1.0));

    gradN(0,3)  = 0.25*(-(u-1.0)*(v+1.0) - (v+1.0)*(u - v + 1.0));
    gradN(1,3)  = 0.25*(+(u-1.0)*(v+1.0) - (u-1.0)*(u - v + 1.0));

    gradN(0,4)  = 0.50*(2*u*(v-1.0));
    gradN(1,4)  = 0.50*(u*u-1.0);

    gradN(0,5)  = 0.50*(-v*v + 1.0);
    gradN(1,5)  = 0.50*(-2.0*(u+1.0)*v);

    gradN(0,6)  = 0.50*(-2.0*u*(v+1.0));
    gradN(1,6)  = 0.50*(-u*u + 1.0);

    gradN(0,7)  = 0.50*(v*v - 1.0);
    gradN(1,7)  = 0.50*(2.0*(u-1)*v);

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JQuadElement ::  getXYCoords( const vector<Point2D> &xyCoords, const Point2D &uv, Point2D &xy)
{
    double u = uv[0];
    double v = uv[1];
    double eps = 1.0E-06;

    if( u < -1.0 - eps ||  u > 1.0 + eps) return 1;
    if( v < -1.0 - eps ||  v > 1.0 + eps) return 2;

    if( xyCoords.size() < 4) {
        cout << "Warning: Quad element must have minimum four nodes " << endl;
        return 3;
    }

    vector<double>  N;
    getShapeFunc( uv, N);

    double xsum = 0.0;
    double ysum = 0.0;

    int nSize = N.size();
    for( int i = 0; i < nSize; i++) {
        xsum += xyCoords[i][0]*N[i];
        ysum += xyCoords[i][1]*N[i];
    }

    xy[0] = xsum;
    xy[1] = ysum;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

double JQuadElement :: getJacobian( const vector<Point2D> &xyCoords, const Point2D &uv)
{
    int npoints = xyCoords.size();
    MatrixXd gradN, xy;

    xy.resize(npoints,2);
    for( int i = 0; i < npoints; i++) {
        xy.coeffRef(i,0) = xyCoords[i][0];
        xy.coeffRef(i,1) = xyCoords[i][1];
    }

    getShapeDeriv(uv, gradN);
    MatrixXd J  = gradN*xy;
    double dJ =  J(0,0)*J(1,1) - J(0,1)*J(1,0);
    return dJ;
}

///////////////////////////////////////////////////////////////////////////////

double JQuadElement :: eval( const Point2D &uv, const vector<double> &phi)
{
    vector<double> N;
    int err = getShapeFunc(uv, N);
    if( err ) return 0.0;

    assert( N.size() == phi.size() );

    int numnodes = phi.size();
    double sum = 0.0;
    for( int i = 0; i < numnodes; i++)
        sum += N[i]*phi[i];
    return sum;
}

/////////////////////////////////////////////////////////////////////////////////
Array2D JQuadElement :: eval( const Point2D &uv, const vector<Array2D> &phi)
{
    Array2D newval;
    vector<double> N;
    int err = getShapeFunc(uv, N);
    if( err ) return newval;

    assert( N.size() == phi.size() );

    int numnodes = phi.size();
    double xsum = 0.0;
    double ysum = 0.0;
    for( int i = 0; i < numnodes; i++) {
        xsum += N[i]*phi[i][0];
        ysum += N[i]*phi[i][0];
    }
    newval[0] = xsum;
    newval[1] = ysum;
    return newval;

}

/////////////////////////////////////////////////////////////////////////////////

Array3D JQuadElement :: eval( const Point2D &uv, const vector<Array3D> &phi)
{
    Array3D newval;
    vector<double> N;
    int err = getShapeFunc(uv, N);
    if( err ) return newval;

    assert( N.size() == phi.size() );

    int numnodes = phi.size();
    double xsum = 0.0;
    double ysum = 0.0;
    double zsum = 0.0;
    for( int i = 0; i < numnodes; i++) {
        xsum += N[i]*phi[i][0];
        ysum += N[i]*phi[i][1];
        zsum += N[i]*phi[i][2];
    }
    newval[0] = xsum;
    newval[1] = ysum;
    newval[2] = zsum;
    return newval;

}
/////////////////////////////////////////////////////////////////////////////////

JFEMSpace :: JFEMSpace()
{
    setOrder(1);
}

///////////////////////////////////////////////////////////////////////////////
void JFEMSpace :: setOrder( int n)
{
    assert(n >= 1 && n <= 2);
    edgeElement.setOrder(n);
    triElement.setOrder(n);
    quadElement.setOrder(n);
}
///////////////////////////////////////////////////////////////////////////////

#ifdef CSV

double FEMSpace :: getQuad1Jacobian( const vector<Point2D> &xyCoords, const Point2D &uv)
{
    MatrixXd gradN, xy;
    xy.resize(4,2);
    for( int i = 0; i < 4; i++) {
        xy.coeffRef(i,0) = xyCoords[i][0];
        xy.coeffRef(i,1) = xyCoords[i][1];
    }

    getQuad1ShapeDeriv(uv, gradN);
    MatrixXd J  = gradN*xy;
    double dJ =  J(0,0)*J(1,1) - J(0,1)*J(1,0);
    return dJ;
}

//////////////////////////////////////////////////////////////////////////////////

void FEMSpace :: getFaceGaussXY(const vector<Point2D> &nodeCoords, vector<Point2D> &xygauss)
{
    xygauss.clear();
    int np = nodeCoords.size();

    int ngauss = 0;
    if( np == 3) ngauss = getNumFaceGaussPoints(TRIANGLE);
    if( np == 4) ngauss = getNumFaceGaussPoints(QUAD);
    if( ngauss == 0) return;

    /*
        xygauss.resize(ngauss);

        vector<double> phi(np);
        // Calculate the "X" Coordinates ..
        for( int i = 0; i < np; i++) phi[i] = nodeCoords[i][0];

        for( int i = 0; i < numFaceGaussPnts; i++) {
            xygauss[i][0] = evalFunc( faceGaussPnts[i][0], faceGaussPnts[i][1], phi);
        }


        // Calculate the "Y" Coordinates ..
        for( int i = 0; i < nsize; i++) phi[i] = nodeCoords[i][1];
        for( int i = 0; i < numFaceGaussPnts; i++)
            xygauss[i][1] = evalFunc( faceGaussPnts[i][0], faceGaussPnts[i][1], phi);
    */

}
///////////////////////////////////////////////////////////////////////////////
/*
void  FEMSpace :: getXYCoords(int elem, const vector<Point2D> &polyPoints, const Point2D &uv, Point2D &xy)
{
    int n = polyPoints.size();

    vector<double>  phi(n);
    for( int i = 0; i < n; i++)
        phi[i] = polyPoints[i][0];
    xy[0] = evalFunc( elem, uv, phi);

    for( int i = 0; i < n; i++)
        phi[i] = polyPoints[i][1];
    xy[1] = evalFunc( elem, uv, phi);
}
*/
///////////////////////////////////////////////////////////////////////////////
void JFEMSpace::getUVCoords(vector<Point2D> &polyPoints, const Point2D &xy, Point2D &uv)
{
    TriGeometry::getUVCoords( &polyPoints[0][0], &polyPoints[1][0], &polyPoints[2][0], &xy[0], &uv[0]);
}
///////////////////////////////////////////////////////////////////////////////

int SelfTangleRegion( const JFacePtr &face, vector<Point2D> &uvCoords, vector<Point2D> &xyCoords)
{
    uvCoords.clear();
    int nx = 100.0;
    int ny = 100.0;

    double dx = 2.0/(double)(nx-1);
    double dy = 2.0/(double)(ny-1);

    Matrix2d  J;
    Point2D uv;
    FEMSpace  femSpace;
    MatrixXd gradN, xy;
//  getCoords(face, xy);

    double dJ;
    for( int j = 0; j < ny; j++) {
        uv[1] = -1.0 + j*dy;
        for( int i = 0; i < nx; i++) {
            uv[0] = -1.0 + i*dx;
            /*
                           femSpace.getShapeDeriv(face, uv, gradN);
                           J  = gradN*xy;
                           dJ =  J(0,0)*J(1,1) - J(0,1)*J(1,0);
                           if( dJ < 0.0) {
                               uvPoints.push_back(uv);
                               xy = femSpace->evalFunc(FEMSpace::QUAD, uv, points);
                           }
            */
        }
    }
}

//////////////////////////////////////////////////////////////////////////////
int JConcaveQuadElement :: getShapeFunc1( int corner, const Point2D &ab, const Point2D &uv, vector<double> &N)
{
    assert( ab[0] >= -1.0 && ab[0] <= 1.0);
    assert( ab[1] >= -1.0 && ab[1] <= 1.0);

    int err;

    switch( corner )
    {
    case 0:
        err = getCorner0Func(ab, uv, N);
        break;
    case 1:
        err = getCorner1Func(ab, uv, N);
        break;
    case 2:
        err = getCorner2Func(ab, uv, N);
        break;
    case 3:
        err = getCorner3Func(ab, uv, N);
        break;
    }
    return err;

}
//////////////////////////////////////////////////////////////////////////////
int JConcaveQuadElement :: getCorner0Func( const Point2D &ab, const Point2D &uv, vector<double> &N)
{
    N.resize(8);
    double a = ab[0];
    double b = ab[1];
    double x = uv[0];
    double y = uv[1];

    N[0] =  (x-1)*(y-1)*(1+x+ y + a*y + b*x - a*b)/((a*a-1)*(b*b-1));
    N[1] = -(a-x)*(b-y)*(y-1)/(2.0*(a-1)*(b+1));
    N[2] =  (a-x)*(b-y)*(x+y)/(2.0*(a-1)*(b-1));
    N[3] = -(a-x)*(x-1)*(b-y)/(2.0*(a+1)*(b-1));

    N[4] = -(x-1)*(b-y)*(y-1)/(2.0*(a-1)*(b+1));
    N[5] = (a-x)*(y-1)*(1 -b + b*x + x + 2*y)/(2.0*(a-1)*(b*b-1));
    N[6] =  (x-1)*(b-y)*(1+2*x -a + a*y + y)/(2.0*(a*a-1)*(b-1));
    N[7] = -(a-x)*(x-1)*(y-1)/(2.0*(a+1)*(b-1));
    return 0;
}
//////////////////////////////////////////////////////////////////////////////

int JConcaveQuadElement :: getCorner1Func( const Point2D &ab, const Point2D &uv, vector<double> &N)
{
    N.resize(8);
    double a = ab[0];
    double b = ab[1];
    double x = uv[0];
    double y = uv[1];

    N[0] = -(a-x)*(b-y)*(y-1)/(2.0*(a+1)*(b+1));
    N[1] =  (1+x)*(y-1)*(-1+x + b*x - y - a*b + a*y)/((a-1)*(a+1)*(b-1)*(b+1));
    N[2] =  (a-x)*(x+1)*(b-y)/(2.0*(a-1)*(b-1));
    N[3] =  (a-x)*(x-y)*(y-b)/(2.0*(a+1)*(b-1));

    N[4] = -(x+1)*(b-y)*(y-1)/(2.0*(a+1)*(b+1));
    N[5] = (a-x)*(1+x)*(y-1)/( 2.0*(a-1)*(b-1));
    N[6] = (x+1)*(b-y)*(-1 + 2*x -a + a*y -y)/(2.0*(a-1)*(a+1)*(b-1));
    N[7] = -(a-x)*(x-1)*(y-1)/(2.0*(a+1)*(b-1));
    return 0;
}

//////////////////////////////////////////////////////////////////////////////

int JConcaveQuadElement :: getCorner2Func( const Point2D &ab, const Point2D &uv, vector<double> &N)
{
    N.resize(8);
    double a = ab[0];
    double b = ab[1];
    double x = uv[0];
    double y = uv[1];

    N[0] = -(a-x)*(b-y)*(x+y)/(2.0*(a+1)*(b+1));
    N[1] =  (a-x)*(x+1)*(b-y)/(2.0*(a-1)*(b+1));
    N[2] = -(x+1)*(y+1)*(-1+x +y- b*x -a*y + a*b)/((a*a-1)*(b*b-1));
    N[3] =  (a-x)*(b-y)*(y+1)/(2.0*(a+1)*(b-1));

    N[4] =  (1+x)*(b-y)*(-1+ 2*x + y -a - a*y)/(2.0*(a*a-1)*(b+1));
    N[5] =  (a-x)*(1+x)*(1+y)/(2.0*(a-1)*(b+1));
    N[6] =  (1+x)*(b-y)*(y+1)/(2.0*(a+1)*(b-1));
    N[7] = -(a-x)*(1+b-x+b*x -2*y)*(y+1)/(2.0*(a+1)*(b*b-1));

    return 0;
}
//////////////////////////////////////////////////////////////////////////////

int JConcaveQuadElement :: getCorner3Func( const Point2D &ab, const Point2D &uv, vector<double> &N)
{
    N.resize(8);
    double a = ab[0];
    double b = ab[1];
    double x = uv[0];
    double y = uv[1];

    N[0] = -(a-x)*(x-1)*(b-y)/(2.0*(a+1)*(b+1));
    N[1] =  (a-x)*(b-y)*(x-y)/(2.0*(a-1)*(b+1));
    N[2] =  (a-x)*(b-y)*(y+1)/(2.0*(a-1)*(b-1));
    N[3] =  (x-1)*(y+1)*(-1.0 -x +b*x + y -a*b  + a*y)/((a*a-1)*(b*b-1));

    N[4] =  (x-1)*(b-y)*(-1.0 + a - 2.*x + y + a*y)/(2.0*(a*a-1)*(b+1));
    N[5] =  (a-x)*(1+y)*(-1-b + b*x - x + 2*y)/(2.0*(a-1)*(1+b*b));
    N[6] =  (x-1)*(b-y)*(y+1)/(2.0*(a-1)*(b-1));
    N[7] = -(a-x)*(x-1)*(y+1)/(2.0*(a+1)*(b-1));
    return 0;
}

///////////////////////////////////////////////////////////////////////////////
int JConcaveQuadElement :: getConcaveCorner( const vector<Point2D> &xyPoints)
{
    double J;
    double eps = -1.0E-10;

    Point2D uv;

    uv[0] = -1.0;
    uv[1] = -1.0;
    J = stdQuad.getJacobian( xyPoints, uv );
    if( J < eps) return 0;

    uv[0] =  1.0;
    uv[1] = -1.0;
    J = stdQuad.getJacobian( xyPoints, uv);
    if( J < eps) return 1;

    uv[0] =  1.0;
    uv[1] =  1.0;
    J = stdQuad.getJacobian( xyPoints, uv);
    if( J < eps) return 2;

    uv[0] = -1.0;
    uv[1] =  1.0;
    J = stdQuad.getJacobian( xyPoints, uv);
    if( J < eps) return 3;

    return -1;
}
///////////////////////////////////////////////////////////////////////////////
void JConcaveQuadElement :: zeroSearch( const vector<Point2D> &xyPoints, Point2D &p0, Point2D &p1, Point2D &pm, int &found)
{
    if( found ) return;

    double J0 = stdQuad.getJacobian( xyPoints, p0 );
    if( fabs(J0) < 1.0E-10) {
        pm = p0;
        found = 1;
        return;
    }

    double J1 = stdQuad.getJacobian( xyPoints, p1 );
    if( fabs(J1) < 1.0E-10) {
        pm = p1;
        found = 1;
        return;
    }

    if( J0*J1 > 0.0) return;

    pm[0] = 0.5*( p0[0] + p1[0] );
    pm[1] = 0.5*( p0[1] + p1[1] );

    double Jm = stdQuad.getJacobian( xyPoints, pm );
    if( fabs(Jm) < 1.0E-10) {
        found = 1;
        return;
    }

    if( J0*Jm < 0.0) {
        p1 = pm;
        zeroSearch(xyPoints, p0, p1, pm, found);
        return;
    }


    if( J1*Jm < 0.0) {
        p0 = pm;
        zeroSearch(xyPoints, p0, p1, pm, found);
        return;
    }
}
///////////////////////////////////////////////////////////////////////////////

int JConcaveQuadElement :: getZeroJacobianLine( const vector<Point2D> &xyPoints, Point2D &pfirst, Point2D &psecond)
{
    Point2D p0, p1, pm;
    Point2D points[2];
    int err, found, index = 0;

    p0[0] =  -1.0;
    p0[1] =  -1.0;
    p1[0] =   1.0;
    p1[1] =  -1.0;
    found = 0;
    zeroSearch( xyPoints, p0, p1, pm, found);
    if( found )  {
        points[index] = pm;
        index++;
    }

    p0[0] =   1.0;
    p0[1] =  -1.0;
    p1[0] =   1.0;
    p1[1] =   1.0;
    found = 0;
    zeroSearch( xyPoints, p0, p1, pm, found);
    if( found )  {
        points[index] = pm;
        index++;
    }
    if( index == 2 ) {
        pfirst  = points[0];
        psecond = points[1];
        return 0;
    }

    p0[0] =   1.0;
    p0[1] =   1.0;
    p1[0] =  -1.0;
    p1[1] =   1.0;
    found = 0;
    zeroSearch( xyPoints, p0, p1, pm, found);
    if( found )  {
        points[index] = pm;
        index++;
    }
    if( index == 2 ) {
        pfirst  = points[0];
        psecond = points[1];
        return 0;
    }

    p0[0] =  -1.0;
    p0[1] =   1.0;
    p1[0] =  -1.0;
    p1[1] =  -1.0;
    found = 0;
    zeroSearch( xyPoints, p0, p1, pm, found);
    if( found )  {
        points[index] = pm;
        index++;
    }
    if( index == 2 ) {
        pfirst  = points[0];
        psecond = points[1];
        return 0;
    }

    return 1;
}
///////////////////////////////////////////////////////////////////////////////
int JConcaveQuadElement :: getNegativeJacobianRegion( const vector<Point2D> &xyPoints, vector<Point2D> &uvCoords)
{
    /*
        uvCoords.clear();
        xyCoords.clear();

        int nx = 100.0;
        int ny = 100.0;

        double dx = 2.0/(double)(nx-1);
        double dy = 2.0/(double)(ny-1);

        Matrix2d  J;
        Point2D uv;
        MatrixXd gradN, xy;

        getCoords(face, xy);
        Point2D p2d;

        double dJ;
        for( int j = 0; j < ny; j++) {
            uv[1] = -1.0 + j*dy;
            for( int i = 0; i < nx; i++) {
                uv[0] = -1.0 + i*dx;
                getShapeDeriv(face, uv, gradN);
                J  = gradN*xy;
                dJ =  J(0,0)*J(1,1) - J(0,1)*J(1,0);
                if( dJ <= 0.0) {
                  getXYCoords(face, uv, p2d);
                    uvCoords.push_back(uv);
                    xyCoords.push_back(p2d);
                }
            }
        }
    */
}

///////////////////////////////////////////////////////////////////////////////
int JConcaveQuadElement :: getReflection0( const vector<Point2D> &xyPoints, Point2D &ab)
{
    if( xyPoints.size() != 4 ) {
        cout << "Warning : Self tangling happens only non-simplicial element" << endl;
        return 1;
    }

    Point2D p0, p1;
    int found;

    p0[0] =  -1.0;
    p0[1] =  -1.0;
    p1[0] =   1.0;
    p1[1] =  -1.0;
    found = 0;
    Point2D p4;
    zeroSearch( xyPoints, p0, p1, p4, found);
    if( !found ) return 2;

    Point2D p7;
    p0[0] =  -1.0;
    p0[1] =  -1.0;
    p1[0] =  -1.0;
    p1[1] =   1.0;
    found = 0;
    zeroSearch( xyPoints, p0, p1, p7, found);
    if( !found ) return 2;

    ab[0] = p4[0];
    ab[1] = p7[1];

#ifdef DEBUG
    QuadElement stdQuad;
    stdQuad.getXYCoord(xyPoints[0], ab, xy);
    double dist = JJMath::length(xyPoints, xy);
    if( dist > 1.E-10)
        cout << "Warning: Reflection point do not match " << dist << endl;
#endif

    return 0;
}
///////////////////////////////////////////////////////////////////////////////

int JConcaveQuadElement :: splitCorner0( const vector<Point2D> &xyPoints, vector<Point2D> &uvCoords, vector<Array4I> &quads)
{
    Point2D ab;
    getReflection0(xyPoints, ab);

    double a = ab[0];
    double b = ab[1];

    uvCoords.resize(8);

    uvCoords[0][0] =  a;
    uvCoords[0][1] =  b;

    uvCoords[1][0] =  1.0;
    uvCoords[1][1] = -1.0;

    uvCoords[2][0] =  1.0;
    uvCoords[2][1] =  1.0;

    uvCoords[3][0] = -1.0;
    uvCoords[3][1] =  1.0;

    uvCoords[4][0] = a;
    uvCoords[4][1] = -1.0;

    uvCoords[5][0] =  1.0;
    uvCoords[5][1] =  b;

    uvCoords[6][0] =  a;
    uvCoords[6][1] =  1.0;

    uvCoords[7][0] = -1.0;
    uvCoords[7][1] =  b;

    quads.resize(3);
    quads[0][0] = 1;
    quads[0][1] = 5;
    quads[0][2] = 0;
    quads[0][3] = 4;

    quads[1][0] = 2;
    quads[1][1] = 6;
    quads[1][2] = 0;
    quads[1][3] = 5;

    quads[2][0] = 3;
    quads[2][1] = 7;
    quads[2][2] = 0;
    quads[2][3] = 6;

    return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////

int JConcaveQuadElement :: getReflection1( const vector<Point2D> &xyPoints, Point2D &ab)
{
    if( xyPoints.size() != 4 ) {
        cout << "Warning : Self tangling happens only non-simplicial element" << endl;
        return 1;
    }

    Point2D p0, p1;
    int found;

    p0[0] =  -1.0;
    p0[1] =  -1.0;
    p1[0] =   1.0;
    p1[1] =  -1.0;
    found = 0;
    Point2D p4;
    zeroSearch( xyPoints, p0, p1, p4, found);
    if( !found ) return 2;

    p0[0] =   1.0;
    p0[1] =  -1.0;
    p1[0] =   1.0;
    p1[1] =   1.0;
    found = 0;
    Point2D p5;
    zeroSearch( xyPoints, p0, p1, p5, found);
    if( !found ) return 2;

    ab[0] = p4[0];
    ab[1] = p5[1];

    return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////

int JConcaveQuadElement :: splitCorner1( const vector<Point2D> &xyPoints, vector<Point2D> &uvCoords, vector<Array4I> &quads)
{
    Point2D ab;
    getReflection1( xyPoints, ab);

    double a = ab[0];
    double b = ab[1];

    uvCoords.resize(8);
    uvCoords[0][0] = -1.0;
    uvCoords[0][1] = -1.0;

    uvCoords[1][0] =  a;
    uvCoords[1][1] =  b;

    uvCoords[2][0] =  1.0;
    uvCoords[2][1] =  1.0;

    uvCoords[3][0] = -1.0;
    uvCoords[3][1] =  1.0;

    uvCoords[4][0] = a;
    uvCoords[4][1] = -1.0;

    uvCoords[5][0] = 1.0;
    uvCoords[5][1] = b;

    uvCoords[6][0] = a;
    uvCoords[6][1] =  1.0;

    uvCoords[7][0] = -1.0;
    uvCoords[7][1] =  b;

    quads.resize(3);
    quads[0][0] = 0;
    quads[0][1] = 4;
    quads[0][2] = 1;
    quads[0][3] = 7;

    quads[1][0] = 2;
    quads[1][1] = 6;
    quads[1][2] = 1;
    quads[1][3] = 5;

    quads[2][0] = 3;
    quads[2][1] = 7;
    quads[2][2] = 1;
    quads[2][3] = 6;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////

int JConcaveQuadElement :: getReflection2( const vector<Point2D> &xyPoints, Point2D &ab)
{
    if( xyPoints.size() != 4 ) {
        cout << "Warning : Self tangling happens only non-simplicial element" << endl;
        return 1;
    }

    Point2D p0, p1;
    int found;

    p0[0] =   1.0;
    p0[1] =  -1.0;
    p1[0] =   1.0;
    p1[1] =   1.0;
    found = 0;
    Point2D p5;
    zeroSearch( xyPoints, p0, p1, p5, found);
    if( !found ) return 2;

    p0[0] =  1.0;
    p0[1] =  1.0;
    p1[0] =  -1.0;
    p1[1] =   1.0;
    found = 0;
    Point2D p6;
    zeroSearch( xyPoints, p0, p1, p6, found);
    if( !found ) return 2;

    ab[0] = p6[0];
    ab[1] = p5[1];
    return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////

int JConcaveQuadElement :: splitCorner2( const vector<Point2D> &xyPoints, vector<Point2D> &uvCoords, vector<Array4I> &quads)
{
    Point2D ab;
    getReflection2( xyPoints, ab);
    double a = ab[0];
    double b = ab[1];

    uvCoords.resize(8);

    uvCoords[0][0] = -1.0;
    uvCoords[0][1] = -1.0;

    uvCoords[1][0] =  1.0;
    uvCoords[1][1] = -1.0;

    uvCoords[2][0] =  a;
    uvCoords[2][1] =  b;

    uvCoords[3][0] = -1.0;
    uvCoords[3][1] =  1.0;

    uvCoords[4][0] =  a;
    uvCoords[4][1] = -1.0;

    uvCoords[5][0] =  1.0;
    uvCoords[5][1] =  b;

    uvCoords[6][0] =  a;
    uvCoords[6][1] =  1.0;

    uvCoords[7][0] = -1.0;
    uvCoords[7][1] =  b;

    quads.resize(3);
    quads[0][0] = 0;
    quads[0][1] = 4;
    quads[0][2] = 2;
    quads[0][3] = 7;

    quads[1][0] = 1;
    quads[1][1] = 5;
    quads[1][2] = 2;
    quads[1][3] = 4;

    quads[2][0] = 3;
    quads[2][1] = 7;
    quads[2][2] = 2;
    quads[2][3] = 6;

    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////

int JConcaveQuadElement :: getReflection3( const vector<Point2D>  &xyPoints, Point2D &ab)
{
    if( xyPoints.size() != 4 ) {
        cout << "Warning : Self tangling happens only non-simplicial element" << endl;
        return 1;
    }

    Point2D p0, p1;
    int found;

    p0[0] =   1.0;
    p0[1] =   1.0;
    p1[0] =  -1.0;
    p1[1] =   1.0;
    found = 0;
    Point2D p6;
    zeroSearch( xyPoints, p0, p1, p6, found);
    if( !found ) return 2;

    p0[0] =  -1.0;
    p0[1] =   1.0;
    p1[0] =  -1.0;
    p1[1] =  -1.0;
    found = 0;
    Point2D p7;
    zeroSearch( xyPoints, p0, p1, p7, found);
    if( !found ) return 2;

    ab[0] = p6[0];
    ab[1] = p7[1];
    return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////

int JConcaveQuadElement :: splitCorner3( const vector<Point2D>  &xyPoints, vector<Point2D> &uvCoords, vector<Array4I> &quads)
{
    Point2D ab;
    getReflection3( xyPoints, ab);

    double a = ab[0];
    double b = ab[1];

    uvCoords.resize(8);
    uvCoords[0][0] = -1.0;
    uvCoords[0][1] = -1.0;

    uvCoords[1][0] =  1.0;
    uvCoords[1][1] = -1.0;

    uvCoords[2][0] =  1.0;
    uvCoords[2][1] =  1.0;

    uvCoords[3][0] =  a;
    uvCoords[3][1] =  b;

    uvCoords[4][0] =  a;
    uvCoords[4][1] = -1.0;

    uvCoords[5][0] = 1.0 ;
    uvCoords[5][1] = b;

    uvCoords[6][0] =  a;
    uvCoords[6][1] =  1.0;

    uvCoords[7][0] = -1.0;
    uvCoords[7][1] =  b;

    quads.resize(3);
    quads[0][0] = 0;
    quads[0][1] = 4;
    quads[0][2] = 3;
    quads[0][3] = 7;

    quads[1][0] = 1;
    quads[1][1] = 5;
    quads[1][2] = 3;
    quads[1][3] = 4;

    quads[2][0] = 2;
    quads[2][1] = 6;
    quads[2][2] = 3;
    quads[2][3] = 5;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////

int JConcaveQuadElement :: split( const vector<Point2D> &xyPoints, vector<Point2D> &uvCoords, vector<Array4I> &quads)
{
    Point2D uv;

    uv[0] = -1.0;
    uv[1] = -1.0;
    double J0 = stdQuad.getJacobian(xyPoints, uv);
    if( J0 < 0.0) return splitCorner0( xyPoints, uvCoords, quads);

    uv[0] = 1.0;
    uv[1] = -1.0;
    double J1 = stdQuad.getJacobian(xyPoints, uv);
    if( J1 < 0.0) return splitCorner1( xyPoints, uvCoords, quads);

    uv[0] = 1.0;
    uv[1] = 1.0;
    double J2 = stdQuad.getJacobian(xyPoints, uv);
    if( J2 < 0.0) return splitCorner2( xyPoints, uvCoords, quads);

    uv[0] =-1.0;
    uv[1] = 1.0;
    double J3 = stdQuad.getJacobian(xyPoints, uv);
    if( J3 < 0.0) return splitCorner3( xyPoints, uvCoords, quads);
    return 1;

}
///////////////////////////////////////////////////////////////////////////////////////////////

int JConcaveQuadElement :: getReflectionPoint( const vector<Point2D> &xyPoints, Point2D &ab)
{
    int c = getConcaveCorner( xyPoints );
    if( c < 0) return 1;
    switch(c)
    {
    case 0:
        return getReflection0(xyPoints,ab);
    case 1:
        return getReflection1(xyPoints,ab);
    case 2:
        return getReflection2(xyPoints,ab);
    case 3:
        return getReflection3(xyPoints,ab);
    }
    return 1;
}

///////////////////////////////////////////////////////////////////////////////////////////////

int JConcaveQuadElement :: getCorner0Deriv(const Point2D &ab, const Point2D &uv, MatrixXd  &gradN)
{
    double a = ab[0];
    double b = ab[1];

    double x = uv[0];
    double y = uv[1];

    gradN.resize(2,4);

    gradN(0,0) = (y-1)*( 2*x + y + a*y - b*(1+a-2*x))/((a*a-1)*(b*b-1));
    gradN(1,0) = (x-1)*( x + b*x + 2*y -a*(1+b-2*y))/((a*a-1)*(b*b-1));

    gradN(0,1) = (b-y)*(y-1)/(2.0*(a-1)*(b-1));
    gradN(1,1) = (x-a)*(1+b-2*y)/(2.0*(a-1)*(b+1));

    gradN(0,2) = (b-y)*(a-2*x-y)/(2.0*(a-1)*(b-1));
    gradN(0,3) = (a-x)*(b-x-2*y)/(2.0*(a-1)*(b-1));

    gradN(1,2) = (y-b)*(1+a-2*x)/(2.0*(a+1)*(b-1));
    gradN(1,3) = (a-x)*(x-1)/(2.0*(a+1)*(b-1));

    return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////

int JConcaveQuadElement :: getCorner1Deriv(const Point2D &ab, const Point2D &uv, MatrixXd  &gradN)
{
    double a = ab[0];
    double b = ab[1];

    double x = uv[0];
    double y = uv[1];

    gradN.resize(2,4);

    gradN(0,0) = (b-y)*(y-1)/(2.0*(a+1)*(b+1));
    gradN(1,0) = (x-a)*(1+b-2*y)/(2.0*(a+1)*(b+1));

    gradN(0,1) = (y-1)*(b-a*b + 2*x + 2*b*x - y + a*y)/((a*a-1)*(b*b-1));
    gradN(1,1) = (x+1)*((b+1)*(x-a) + 2*y*(a-1))/((a*a-1)*(b*b-1));

    gradN(0,2) = (b-y)*(-1 + a - 2*x)/(2.0*(a-1)*(b-1));
    gradN(0,3) = (x-a)*(x+1)/(2.0*(a-1)*(b-1));

    gradN(1,2) = (y-b)*(a-2*x+y)/(2.0*(a+1)*(b-1));
    gradN(1,3) = (a-x)*(b+x-2*y)/(2.0*(a+1)*(b-1));

    return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////

int JConcaveQuadElement :: getCorner2Deriv(const Point2D &ab, const Point2D &uv, MatrixXd  &gradN)
{
    double a = ab[0];
    double b = ab[1];

    double x = uv[0];
    double y = uv[1];

    gradN.resize(2,4);

    gradN(0,0) = (b-y)*(-a+2*x + y)/(2*(a+1)*(1+b));
    gradN(1,0) = (a-x)*(-b+x+2*y)/( 2*(a+1)*(b+1));

    gradN(0,1) = (b-y)*(-1+a-2*x)/(2*(a-1)*(b+1));
    gradN(1,1) = (x-a)*(x+1)/(2*(a-1)*(b+1));

    gradN(0,2) = (y+1)*(b -a*b - 2*x + 2*b*x -y + a*y)/((a*a-1)*(b*b-1));
    gradN(0,3) = (x+1)*((b-1)*(x-a) -2*y + 2*a*y)/((a*a-1)*(b*b-1));

    gradN(1,2) = (y-b)*(y+1)/(2*(a+1)*(b-1));
    gradN(1,3) = (a-x)*(-1+b-2*y)/(2*(a+1)*(b-1));

    return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////

int JConcaveQuadElement :: getCorner3Deriv(const Point2D &ab, const Point2D &uv, MatrixXd  &gradN)
{
    double a = ab[0];
    double b = ab[1];

    double x = uv[0];
    double y = uv[1];

    gradN.resize(2,4);

    gradN(0,0) = (y-b)*(1+a-2*x)/(2*(a+1)*(b+1));
    gradN(1,0) = (a-x)*(x-1)/(2*(a+1)*(b+1));

    gradN(0,1) = (b-y)*(a-2*x+y)/(2*(a-1)*(b+1));
    gradN(1,1) = (x-a)*(b+x-2*y)/(2*(a-1)*(b+1));

    gradN(0,2) = (y-b)*(y+1)/(2*(a-1)*(b-1));
    gradN(0,3) = (a-x)*(-1 + b -2*y)/(2*(a-1)*(b-1));

    gradN(1,2) = (y+1)*(-b*(1+a-2*x)-2*x+y+a*y)/((a*a-1)*(b*b-1));
    gradN(1,3) = (x-1)*((1-b)*(a-x) + 2*(a+1)*y)/((a*a-1)*(b*b-1));

    return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////
int JConcaveQuadElement :: getGaussUVCoords( const vector<Point2D> &xyCoords, vector<Point3D> &gaussPoints, int n)
{
    gaussPoints.clear();

    int corner = getConcaveCorner(xyCoords);
    if( corner < 0) return 1;

    vector<Point2D> uvCoords;
    vector<Array4I> quads;
    int err = split( xyCoords, uvCoords, quads);
    if( err ) return err;

    QuadElement  stdQuad;

    vector<Point2D> qPoints(4);
    vector<Point3D> gPoints;

    Point2D uv, xy;
    Point3D uvg;

    double sumarea =  0.0;
    for( int i = 0; i < 3; i++) {
        for( int j = 0; j < 4; j++) {
            int id = quads[i][j];
            qPoints[j] =  uvCoords[id];
        }
        stdQuad.getGaussUVCoords( gPoints);
        double area = JFaceGeometry::getArea2D(qPoints);
        for( int j = 0; j < gPoints.size(); j++) {
            uv[0] = gPoints[j][0];
            uv[1] = gPoints[j][1];
            stdQuad.getXYCoords( qPoints, uv, xy);
            double J = stdQuad.getJacobian(qPoints, uv);
            uvg[0] = xy[0];
            uvg[1] = xy[1];
            uvg[2] = J*gPoints[j][2];   // Scaling ....
            gaussPoints.push_back(uvg);
        }
        sumarea += area;
    }

#ifdef DEBUG
    double wsum = 0.0;
    for(int i = 0; i < 12; i++)
        wsum += gaussPoints[i][2];
    assert( fabs(wsum - sumarea) < 1.0E-06);
#endif

    return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////
#endif
