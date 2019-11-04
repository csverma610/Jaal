#include "FEM2D.hpp"

////////////////////////////////////////////////////////////////////////////////
JFEM2D :: JFEM2D()
{
    shapeOrder  = 1;
    dofPerNode  = 0;
    verbose     = 0;
    numTangledPairs = 0;
    detectTangle     = 1;
    absJacobian  = 0;
    femSpace.reset( new JFEMSpace );
    tangleSpace.reset( new JFEMSpace );
}
////////////////////////////////////////////////////////////////////////////////

void JFEM2D :: setOrder(int order)
{
    if( order < 1) {
        cout << "Info: Setting to minimum order 1 " << endl;
        order = 1;
    }

    if( order > 2) {
        cout << "Info: Setting to maximim order to 2" << endl;
        order = 2;
    }

    shapeOrder = order;
    femSpace->setOrder(order);
}

////////////////////////////////////////////////////////////////////////////////

void JFEM2D :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    globalStatus = initMesh();
}

////////////////////////////////////////////////////////////////////////////////
int JFEM2D :: getNumOfNodes( const JEdgePtr &edge) const
{
    return femSpace->edgeElement.getNumNodes();
}

////////////////////////////////////////////////////////////////////////////////
int JFEM2D :: getNumOfNodes( const JFacePtr &face) const
{
    if( face->getSize(0) == 3 ) return femSpace->triElement.getNumNodes();
    if( face->getSize(0) == 4 ) return femSpace->quadElement.getNumNodes();
    return 0;
}
////////////////////////////////////////////////////////////////////////////////

void JFEM2D :: getNodes( const JFacePtr &face, vector<int> &nodes)
{
    // Give the nodes of elements including "Steiner nodes" on edges ...
    nodes.clear();

    for( int i = 0; i < face->getSize(0); i++)
        nodes.push_back(face->getNodeAt(i)->getID());

    JNodePtr vtx;
    if( shapeOrder == 2) {
        for( int i = 0; i < face->getSize(1); i++) {
            const JEdgePtr &edge = face->getEdgeAt(i);
            edge->getAttribute("Order2Node", vtx);
            nodes.push_back(vtx->getID());
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void JFEM2D :: getShapeFunc( const JFacePtr &face, const Point2D &uv, vector<double> &N)
{
    int nn = face->getSize(0);
    if( nn == 3 ) femSpace->triElement.getShapeFunc( uv, N);
    if( nn == 4 ) femSpace->quadElement.getShapeFunc( uv, N);
}

////////////////////////////////////////////////////////////////////////////////

void JFEM2D :: getShapeDeriv( const JFacePtr &face, const Point2D &uv, MatrixXd &gradN)
{
    int nn = face->getSize(0);
    if( nn == 3 )  femSpace->triElement.getShapeDeriv( uv, gradN);
    if( nn  == 4 ) femSpace->quadElement.getShapeDeriv( uv, gradN);
}

////////////////////////////////////////////////////////////////////////////////
void JFEM2D :: getCoords( const JFacePtr &face, MatrixXd &xy)
{
    int nnodes = face->getSize(0);
    int nrows  = nnodes;
    if( shapeOrder == 2 ) nrows *= 2;

    xy.resize(nrows,2);

    int irow = 0;
    for( int i = 0; i < nnodes; i++) {
        const Point3D &xyz = face->getNodeAt(i)->getXYZCoords();
        xy(irow,0) = xyz[0];
        xy(irow,1) = xyz[1];
        irow++;
    }

    JNodePtr vtx;
    if( shapeOrder == 2) {
        for( int i = 0; i < nnodes; i++) {
            const JEdgePtr &edge = face->getEdgeAt(i);
            edge->getAttribute("Order2Node", vtx);
            const Point3D &xyz = vtx->getXYZCoords();
            xy(irow,0) = xyz[0];
            xy(irow,1) = xyz[1];
            irow++;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
void JFEM2D :: getXYCoords( const JFacePtr &face, const Point2D &uv, Point2D &xy)
{
    /*
        assert( femSpace );
        MatrixXd nodeCoords;
        getCoords(face, nodeCoords);

        int nnodes = nodeCoords.rows();
        vector<double> x(nnodes), y(nnodes);

        for( int i = 0; i < nnodes; i++) {
            x[i] = nodeCoords.coeff(i,0);
            y[i] = nodeCoords.coeff(i,1);
        }

        if( face->getSize(0) == 3 ) {
            xy[0] = femSpace->evalFunc( FEMSpace::TRIANGLE, uv, x);
            xy[1] = femSpace->evalFunc( FEMSpace::TRIANGLE, uv, y);
        }

        if( face->getSize(0) == 4 ) {
            xy[0] = femSpace->evalFunc( FEMSpace::QUAD, uv, x);
            xy[1] = femSpace->evalFunc( FEMSpace::QUAD, uv, y);
            return;
        }
    */
}
////////////////////////////////////////////////////////////////////////////////

void JFEM2D :: mapUV(int triangle, const Point2D &uvtri, Point2D &uvquad)
{
    // Given (u,v) coordinate on a triangle face (in uvtri), Get the (u,v) is the quad mesh.
    // We split the quad into two triangles and then assign the (uq,vq)

    double u[3], v[3];
    switch( triangle)
    {
    case 0:
        u[0] = -1.0;
        u[1] =  1.0;
        u[2] =  1.0;
        v[0] = -1.0;
        v[1] = -1.0;
        v[2] =  1.0;
        break;
    case 1:
        u[0] = -1.0;
        u[1] =  1.0;
        u[2] = -1.0;
        v[0] = -1.0;
        v[1] =  1.0;
        v[2] =  1.0;
        break;
    }

    vector<double> N;
//  femSpace->getShapeFunc(FEMSpace::TRIANGLE, uvtri[0], uvtri[1], N);

    assert(N.size() == 3);

    double usum = 0.0;
    double vsum = 0.0;
    for( int i = 0; i < 3; i++) {
        usum += N[i]*u[i];
        vsum += N[i]*v[i];
    }
    uvquad[0] = usum;
    uvquad[1] = vsum;
    assert( usum >= -1.0 && usum <= 1.0);
    assert( vsum >= -1.0 && vsum <= 1.0);
}

////////////////////////////////////////////////////////////////////////////////

void JFEM2D :: getAssemblyIndices(const JFacePtr &face, vector<int> &globalPos)
{
    //
    // Find the global pos of the nodes on a face. We follow the convention that
    // dof are specified as (u_1, v_1, u_2, v_2 ..... u_n, v_n)
    //
    // Also, first we enumerate the nodes of the elements and then the nodes on
    // the edges ( which are higher order nodes)..
    //
    globalPos.clear();

    int nn = face->getSize(0);
    for( int i = 0; i < nn; i++) {
        int vid = face->getNodeAt(i)->getID();
        for( int j = 0; j < dofPerNode; j++)
            globalPos.push_back(dofPerNode*vid+j);
    }

    JNodePtr vtx;
    if( shapeOrder == 2) {
        for( int i = 0; i < nn; i++) {
            const JEdgePtr &edge = face->getEdgeAt(i);
            edge->getAttribute("Order2Node", vtx);
            int vid = vtx->getID();
            for( int j = 0; j < dofPerNode; j++) {
                globalPos.push_back(dofPerNode*vid+j);
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////

void JFEM2D :: getAssemblyIndices(const JEdgePtr &edge, vector<int> &globalPos)
{
    // Find the global pos of the nodes on a edge. We follow the convention that
    // dof are specified as (u_1, v_1, u_2, v_2 ..... u_n, v_n)
    //
    // Also, first we enumerate the end nodes of a edge  and then any higher
    // order nodes on the edge..

    int vid;
    globalPos.clear();

    vid = edge->getNodeAt(0)->getID();
    for( int j = 0; j < dofPerNode; j++)
        globalPos.push_back(dofPerNode*vid+j);

    vid = edge->getNodeAt(1)->getID();
    for( int j = 0; j < dofPerNode; j++)
        globalPos.push_back(dofPerNode*vid+j);

    JNodePtr vtx;
    if( shapeOrder == 2) {
        edge->getAttribute("Order2Node", vtx);
        vid =  vtx->getID();
        for( int j = 0; j < dofPerNode; j++)
            globalPos.push_back(dofPerNode*vid+j);
    }
}

//////////////////////////////////////////////////////////////////////////////////////

void JFEM2D :: genQuadraticNodes()
{
    if( mesh == nullptr ) return;

    if( verbose )
        cout << "Info: Generating Quadratic mesh ... " << endl;

    // Acutally, I would not this as below, but this is how matlab finds uniques
    // edges and assign new nodes..
    int numNodes = mesh->getSize(0);
    int numFaces = mesh->getSize(2);

    int index = numNodes;
    mesh->deleteEdgeAttribute("Order2Node");
    for( int iface =  0; iface < numFaces; iface++) {
        const JFacePtr &face = mesh->getFaceAt(iface);
        if( face->isActive() ) {
            for( int iedge = 0; iedge < face->getSize(1); iedge++) {
                const JEdgePtr &edge = face->getEdgeAt(iedge);
                if( !edge->hasAttribute("Order2Node") ) {
                    const JNodePtr &v0 = edge->getNodeAt(0);
                    const JNodePtr &v1 = edge->getNodeAt(1);
                    JNodePtr vtx = JNodeGeometry::getMidNode( v0,v1);
                    vtx->setID(index++);
                    edge->setAttribute("Order2Node", vtx);
                }
            }
        }
    }

    if( verbose )
        cout << "Info: Quadratic nodes insertted" << endl;
}

//////////////////////////////////////////////////////////////////////////////////////

int JFEM2D :: initMesh()
{
    if( mesh == nullptr ) return 1;

    if( verbose) {
        stopWatch.reset();
        stopWatch.start();
    }

    if( mesh->getTopology()->getDimension() != 2) {
        cout << "Fatal error: 2D Elasticity is only for triangle and quad mesh " << endl;
        return 1;
    }

    mesh->pruneAll();
    mesh->getTopology()->searchBoundary();

    size_t numBound = mesh->getTopology()->getBoundarySize(1);
    if( numBound == 0) {
        cout << "Fatal error: 2D plane must have boundary edges " << endl;
        return 2;
    }

    if( !mesh->getTopology()->isConsistent() ) {
        cout << "Warning: Initial mesh is not consistent: Making it consistent ..." << endl;
        int err = mesh->getTopology()->getConsistent();
        if( err ) {
            cout << "Fatal error: Could not make the mesh consistent " << endl;
            return 2;
        }
    }

    mesh->enumerate(0);
    mesh->enumerate(1);
    mesh->enumerate(2);

    if( verbose)
        cout << "Detecting negative elements ... " << endl;

    size_t numfaces = mesh->getSize(2);
    for( size_t i = 0;  i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        double area = JFaceGeometry::getSignedArea( face );
        if( fabs(area) < 1.0E-15) {
            cout << "Fatal error: one of the face has zero area " << endl;
            return 1;
        }
    }

    int nCount = 0;
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        short int fsign  = JFaceGeometry::getOrientation2D(face);
        face->setAttribute("Orient", fsign);
        if( fsign < 0.0) nCount++;
    }

    if( nCount == numfaces) {
        cout << "All the faces are negative: Flipping then ... " << endl;
        mesh->getTopology()->reverseAll();
        nCount = 0;
        mesh->deleteFaceAttribute("Orient");
    }

    if( nCount ) {
        cout <<"Info: There are few inverted elements in the mesh ..." << endl;
    }

    nCount = 0;
    for( size_t i = 0;  i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        if( !JFaceGeometry::isSimple(face) ) nCount++;
    }

    if( nCount)
        cout << "Warning : There are some twisted quad elements in the mesh " << endl;

    if( detectTangle ) {
    }

    /*
        if( tangleDetect ) {
            boxFaceID.clear();
            vector<Box> faceBoxes;
            faceBoxes.reserve(numfaces);
            for( size_t i = 0; i < numfaces; i++) {
                const JFacePtr &face = mesh->getFaceAt(i);
                assert( face->getID() == i);
                const BoundingBox  &bx = face->getBoundingBox();
                const Point3D  &pmin   = bx.getLower();
                const Point3D  &pmax   = bx.getUpper();
                faceBoxes.push_back(Bbox(pmin[0], pmin[1], pmax[0], pmax[1]));
                boxFaceID[faceBoxes[i].id()] = i;
            }

            CGAL::box_self_intersection_d( faceBoxes.begin(), faceBoxes.end(),
                                           BoxBoxIntersectCallback);
        }
        if( verbose) {
            stopWatch.stop();
            cout <<"Info: Execution time in detecting tangled boxes " << stopWatch.getSeconds() << endl;
        }
    */

    return 0;
}

//////////////////////////////////////////////////////////////////////////////////////
int  getUVAnalytic(const JFacePtr &face, const Point2D &xy, Point2D &uv)
{
    assert( face->getSize(0) == 4);
    ////////////////////////////////////////////////////////////////////////////////////
    //  Given a Quad element and a specific value, xy. return UV coordinates
    //  using
    //     x = a0 + a1*xi + a2*eta + a3*xi*eta;
    //     y = b0 + b1*xi + b2*eta + b3*xi*eta;
    //
    //     This following is the inteverse of
    //              1    -1   -1   1
    //              1     1   -1  -1
    //              1     1    1   1
    //              1    -1    1  -1
    //
    ////////////////////////////////////////////////////////////////////////////////////
    Matrix4d  mat;
    mat(0,0) = 0.25;
    mat(0,1) = 0.25;
    mat(0,2) = 0.25;
    mat(0,3) = 0.25;
    mat(1,0) = -0.25;
    mat(1,1) =  0.25;
    mat(1,2) = 0.25;
    mat(1,3) = -0.25;
    mat(2,0) = -0.25;
    mat(2,1) = -0.25;
    mat(2,2) = 0.25;
    mat(2,3) =  0.25;
    mat(3,0) =  0.25;
    mat(3,1) = -0.25;
    mat(3,2) = 0.25;
    mat(3,3) = -0.25;

    Vector4d  x, y;
    for( int i = 0; i < 4; i++) {
        const Point3D &xyz = face->getNodeAt(i)->getXYZCoords();
        x(i) = xyz[0];
        y(i) = xyz[1];
    }
    Vector4d a = mat*x;
    Vector4d b = mat*y;

    // Calculate the coefficient of    A*y^2 + By + C = 0;
    //
    double x0 = xy[0] - a[0];
    double y0 = xy[1] - b[0];
    double A  = a[3]*b[2]  - a[2]*b[3];
    double B  = (x0*b[3] + a[1]*b[2])  - (y0*a[3] + a[2]*b[1]);
    double C  = x0*b[1]  - y0*a[1];

    double K = B*B - 4*A*C;
    assert( K >= 0);

    double u = 2.0, v = 2.0;

    if(fabs(A) < 1.0E-15) {
        assert(fabs(B) > 1.0E-10);
        v = -C/B;
    } else {
        double v1 = (-B + sqrt(K))/(2.0*A);
        double v2 = (-B - sqrt(K))/(2.0*A);
        if( v1 >= -1.0 && v1 <= 1.0) v = v1;
        if( v2 >= -1.0 && v2 <= 1.0) v = v2;
    }

    assert( v >= -1.0 && v <= 1.0);

    if( fabs(a[1] + v*a[3]) > 1.0E-10) {
        u = (x0 - a[2]*v)/(a[1] + a[3]*v);
        uv[0] = u;
        uv[1] = v;
        return 0;
    }

    if( fabs(a[3]) > 1.0E-10) {
        v = -a[1]/a[3];
        u = (y0*a[3] + a[1]*b[2])/(a[3]*b[1]-a[1]*b[3]);
        uv[0] =  u;
        uv[1] =  v;
        return 0;
    }

    if( fabs(a[2]) > 1.0E-10) {
        v = x0/a[2];
        u = (y0*a[2] + x0*b[2])/(a[2]*b[1]+x0*b[3]);
        uv[0] =  u;
        uv[1] =  v;
        return 0;
    }

    cout << "Error: Invalid reverse mapping: " << endl;
    exit(0);
    return 1;
}

///////////////////////////////////////////////////////////////////////////////

int  getUVIterative(const JFacePtr &face, const Point2D &xy, Point2D &uv)
{
    assert( face->getSize(0) == 4);
    // Given a Quad face and a xy, return UV values, if the point lies inside it .
    // Using iterative method ..

    Matrix4d  mat;
    mat(0,0) = 0.25;
    mat(0,1) = 0.25;
    mat(0,2) = 0.25;
    mat(0,3) = 0.25;
    mat(1,0) = -0.25;
    mat(1,1) =  0.25;
    mat(1,2) = 0.25;
    mat(1,3) = -0.25;
    mat(2,0) = -0.25;
    mat(2,1) = -0.25;
    mat(2,2) = 0.25;
    mat(2,3) =  0.25;
    mat(3,0) =  0.25;
    mat(3,1) = -0.25;
    mat(3,2) = 0.25;
    mat(3,3) = -0.25;

    Vector4d  x, y;
    for( int i = 0; i < 4; i++) {
        const Point3D &xyz = face->getNodeAt(i)->getXYZCoords();
        x(i) = xyz[0];
        y(i) = xyz[1];
    }
    Vector4d a = mat*x;
    Vector4d b = mat*y;

    Matrix2d  A1;
    A1(0,0) = a[1];
    A1(0,1) = a[2];
    A1(1,0) = b[1];
    A1(1,1) = b[2];
    Matrix2d A = A1.inverse();

    Vector2d param, rhs;
    param[0] = 0.0;
    param[1] = 0.0;
    double xg, yg, dx, dy, dl;
    for( int i = 0; i < 100; i++) {
        rhs[0] =  xy[0] - a[0] - a[3]*param[0]*param[1];
        rhs[1] =  xy[1] - b[0] - b[3]*param[0]*param[1];
        param  =  A*rhs;
        xg = a[0] + a[1]*param[0] + a[2]*param[1] + a[3]*param[0]*param[1];
        yg = b[0] + b[1]*param[0] + b[2]*param[1] + b[3]*param[0]*param[1];
        dx = xy[0] - xg;
        dy = xy[1] - yg;
        dl = sqrt(dx*dx + dy*dy);
        if( dl < 1.0E-12) break;
    }
    if( dl > 1.E0-12)
        cout << "Warning: Iteration did not converge " << dl << endl;

    uv[0] = param[0];
    uv[1] = param[1];
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

bool isSymmetric( const Eigen::SparseMatrix<double> &A)
{
    double eps = 1.0E-10;
    int nrows = A.rows();
    int ncols = A.cols();
    for( int i = 0;   i < nrows; i++) {
        for( int j = i+1; j < ncols; j++) {
            if( fabs(A.coeff(i,j) - A.coeff(j,i)) > eps ) {
                return 0;
            }
        }
    }
    return 1;
}
//
/////////////////////////////////////////////////////////////////////////////////////////
//
bool isSymmetric( const MatrixXd &A)
{
    double eps = 1.0E-10;
    int nrows = A.rows();
    int ncols = A.cols();
    for( int i = 0;   i < nrows; i++) {
        for( int j = i+1; j < ncols; j++) {
            if( fabs(A.coeff(i,j) - A.coeff(j,i)) > eps ) {
                return 0;
            }
        }
    }
    return 1;
}
//
/////////////////////////////////////////////////////////////////////////////////////////
//

void latexMatrix( const MatrixXd &A)
{
    // Print a dense matrix in Latex (useful for the documentation);
    cout << fixed;
    int numRows = A.rows();
    int numCols = A.cols();
    cout << setprecision(5);
    cout << "\\begin{equation}" << endl;
    cout << "\\left[" << endl;
    cout << "\\begin{array}{";
    for( int j = 0; j < numCols; j++)
        cout << "r";
    cout << "}" << endl;
    for( int i = 0; i < numRows; i++) {
        for( int j = 0; j < numCols; j++) {
            cout << A(i,j);
            if( j == numCols-1 )
                cout << "  \\\\ ";
            else
                cout << " & ";
        }
        cout << endl;
    }
    cout << "\\end{array}" << endl;
    cout << "\\right]" << endl;
    cout << "\\end{equation}" << endl;
}
//
///////////////////////////////////////////////////////////////////////////////
//
void latexVector( const vector<double> &v)
{
    //
    // Write a vector in the Latex format (useful for the documentation)
    //
    cout << fixed;
    cout << "\\begin{array}{r}" << endl;
    for( size_t i = 0; i < v.size(); i++)
        cout << v[i] <<  " \\\\" << endl;
    cout << "\\end{array}" << endl;
}
//
///////////////////////////////////////////////////////////////////////////////

void JFEM2D::reduceSystem( Eigen::SparseMatrix<double> &A, vector<double> &b,
                           map<int,double> &fixMap, vector<int> &freeDof)
{
    if( verbose )
        cout << "Info: Reducing the linear system ... " << endl;

    int nrows = A.rows();
    int ncols = A.cols();
    assert( b.size() == size_t(nrows));

    freeDof.resize(nrows);
    for(int i = 0; i < nrows; i++) freeDof[i] = i;

    if( fixMap.empty() ) return;

    for( int i = 0; i < nrows; i++) {
        for( const pair<int,double> &keyVal: fixMap) {
            int j = keyVal.first;
            double val = keyVal.second;
            b[i] = b[i] - A.coeff(i,j)*val;
            A.coeffRef(i,j) = 0.0;
        }
    }

    for( const pair<int,double> &keyVal: fixMap) {
        int i = keyVal.first;
        for( int j = 0; j < ncols; j++) {
            if( A.coeff(i,j) != 0.0) A.coeffRef(i,j) = 0.0;
        }
        A.coeffRef(i,i) = 1.0;
        b[i] = keyVal.second;
    }
    A.prune(0.0);

    /*
        if( debugStep == 5) {
            cout << "Linear System after Dirichlet boundary conditions" << endl;
            latexMatrix(A);
            latexVector(b);
        }
    */

    vector<int> fullDof(nrows);
    for(int i = 0; i < nrows; i++) fullDof[i] = i;

    vector<int> fixDof;
    fixDof.reserve(nrows);
    for( const pair<int,double> &keyVal: fixMap)
        fixDof.push_back(keyVal.first);

    freeDof.clear();
    boost::set_difference(fullDof, fixDof, back_inserter(freeDof));

    vector<int> permrow;
    permrow.reserve(nrows);
    permrow.insert( permrow.end(), freeDof.begin(), freeDof.end() );
    permrow.insert( permrow.end(), fixDof.begin(),  fixDof.end()  );

    Eigen::SparseMatrix<double> permuteMatrix(nrows,nrows);
    for( int i = 0; i < nrows; i++)
        permuteMatrix.coeffRef(i, permrow[i] ) = 1.0;

    A = permuteMatrix*A*permuteMatrix.transpose();

    int nfree = freeDof.size();
    vector<double> bb(nfree);
    for( int i = 0; i < nfree; i++)
        bb[i] = b[permrow[i]];
    b = bb;

    Eigen::SparseMatrix<double> Am(nfree,nfree);
    for( int i = 0; i < nfree; i++) {
        for( int j = 0; j < nfree; j++) {
            double val = A.coeff(i,j);
            if( val != 0.0) Am.insert(i,j) = val;
        }
    }
    A  = Am;

    /*
        if( debugStep == 6) {
            cout << "Reduced linear system after Dirichlet boundary conditions" << endl;
            latexMatrix(A);
            latexVector(b);
        }
    */

    if( verbose )
        cout << "Info: Reducing system ready" << endl;
}

////////////////////////////////////////////////////////////////////////////////

int JFEM2D :: solve()
{
    if( mesh == nullptr ) return 1;
    int topDim = mesh->getTopology()->getDimension();

    if( topDim != 2 ) {
        cout <<"Warning: Poisson 2D is only for triangle and quad elements " << endl;
        return 1;
    }

    numDOF =  mesh->getSize(0);
    if( shapeOrder == 2) {
        cout << "Info: Generting Quadratic nodes " << endl;
        numDOF  +=mesh->getSize(1);
        genQuadraticNodes();
    }

    edgeGaussPoints = femSpace->edgeElement.getGaussUPoints();
    quadGaussPoints = femSpace->quadElement.getGaussUVPoints();
    triGaussPoints  = femSpace->triElement.getGaussUVPoints();

    numNodesPerElement[0] = femSpace->edgeElement.getNumNodes();
    numNodesPerElement[1] = femSpace->triElement.getNumNodes();
    numNodesPerElement[2] = femSpace->quadElement.getNumNodes();

    cout << "#Nodes per edge      " << numNodesPerElement[0] << endl;
    cout << "#Nodes per triangle  " << numNodesPerElement[1] << endl;
    cout << "#Nodes per Quads     " << numNodesPerElement[2] << endl;

    int err = 0;
    err = buildLinearSystem();
    if ( err ) return err;

    err = solveLinearSystem();
    return err;
}

////////////////////////////////////////////////////////////////////////////////
int JFEM2D :: getUVCoords( const JFacePtr &face, const Point2D &xy, Point2D &uv)
{
    int err = 1;
    int numnodes = face->getSize(0);
    cout << "Exit " << endl;
    exit(0);

    /*
        vector<Point2D> points(numnodes);
        for( int i = 0; i < numnodes; i++) {
            const Point3D &xyz = face->getNodeAt(i)->getXYZCoords();
            points[i][0] = xyz[0];
            points[i][1] = xyz[1];
        }

        if( !isInside(face, xy) )  {
            cout << "Fatal error: the point is not inside quad element " << endl;
            exit(0);
        }

        if( numnodes == 3 ) {
            err = TriGeometry::getUVCoords( &points[0][0], &points[1][0], &points[2][0], &xy[0], &uv[0] );
            assert( uv[0] >=0.0 && uv[0] <= 1.0);
            assert( uv[1] >=0.0 && uv[1] <= 1.0);
            return err;
        }

        if( numnodes == 4 ) {
            getTexCoords2(face, xy, uv);
            Point2D xt;
            femshape.getXYCoords( points, uv, xt);
            if( !isInside(face, xt) )  {
                cout << "Fatal error: Recovered point is not inside the quad element " << endl;
                exit(0);
            }
            double dx = xt[0] - xy[0];
            double dy = xt[1] - xy[1];
            double tol = sqrt(dx*dx + dy*dy);
            if( tol > 1.0E-06) {
                cout << "Fatal error: Incorrect (X,Y)->(U,V) mapping in Quad " << endl;
                cout << fixed << setprecision(10) << endl;
                cout << "Quad point difference " << tol << endl;
                cout << points[0][0] << "  " << points[0][1] << endl;
                cout << points[1][0] << "  " << points[1][1] << endl;
                cout << points[2][0] << "  " << points[2][1] << endl;
                cout << points[3][0] << "  " << points[3][1] << endl;
                cout << "Query Point "  << endl;
                cout << xy[0] << "  " << xy[1] << endl;
                cout << "Param Point " << endl;
                cout << uv[0] << "  " << uv[1] << endl;
                cout << "Recover Point " << xt[0] << " " << xt[1] << endl;
                exit(0);
            }
            return 0;
        }
    */
    return err;
}

///////////////////////////////////////////////////////////////////////////////

#ifdef VVDV
double JFEM2D :: getJacobian( const JFacePtr &face, const Point2D &uv)
{
    MatrixXd gradN, xy;
    getCoords(face, xy);

    getShapeDeriv(face, uv, gradN);
    MatrixXd J  = gradN*xy;
    double dJ =  J(0,0)*J(1,1) - J(0,1)*J(1,0);
    return dJ;
}

///////////////////////////////////////////////////////////////////////////////
void JFEM2D :: zeroSearch( const JFacePtr &face, Point2D &p0, Point2D &p1, Point2D &pm, int &found)
{
    if( found ) return;

    double J0 = getJacobian( face, p0 );
    if( fabs(J0) < 1.0E-10) {
        pm = p0;
        found = 1;
        return;
    }
    double J1 = getJacobian( face, p1 );
    if( fabs(J1) < 1.0E-10) {
        pm = p1;
        found = 1;
        return;
    }

    if( J0*J1 > 0.0) return;

    pm[0] = 0.5*( p0[0] + p1[0] );
    pm[1] = 0.5*( p0[1] + p1[1] );

    double Jm = getJacobian( face, pm );
    if( fabs(Jm) < 1.0E-10) {
        found = 1;
        return;
    }

    if( J0*Jm < 0.0) {
        p1 = pm;
        zeroSearch(face, p0, p1, pm, found);
        return;
    }

    if( J1*Jm < 0.0) {
        p0 = pm;
        zeroSearch(face, p0, p1, pm, found);
        return;
    }
}
///////////////////////////////////////////////////////////////////////////////

int JFEM2D :: zeroJacobian( const JFacePtr &face, Point2D &pfirst, Point2D &psecond)
{
    Point2D p0, p1, pm;
    Point2D points[2];
    int err, found, index = 0;

    p0[0] =  -1.0;
    p0[1] =  -1.0;
    p1[0] =   1.0;
    p1[1] =  -1.0;
    found = 0;
    zeroSearch( face, p0, p1, pm, found);
    if( found )  {
        points[index] = pm;
        index++;
    }

    p0[0] =   1.0;
    p0[1] =  -1.0;
    p1[0] =   1.0;
    p1[1] =   1.0;
    found = 0;
    zeroSearch( face, p0, p1, pm, found);
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
    zeroSearch( face, p0, p1, pm, found);
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
    zeroSearch( face, p0, p1, pm, found);
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
int JFEM2D :: getSelfTangleRegion( const JFacePtr &face, vector<Point2D> &uvCoords, vector<Point2D> &xyCoords)
{
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
}

///////////////////////////////////////////////////////////////////////////////
int JFEM2D :: splitSelfTangledCorner0( const JFacePtr &face, vector<Point2D> &uvCoords, vector<Array4I> &quads)
{
    uvCoords.clear();

    if( face->getSize(0) != 4 ) {
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
    zeroSearch( face, p0, p1, p4, found);
    if( !found ) return 2;

    Point2D p7;
    p0[0] =  -1.0;
    p0[1] =  -1.0;
    p1[0] =  -1.0;
    p1[1] =   1.0;
    found = 0;
    zeroSearch( face, p0, p1, p7, found);
    if( !found ) return 2;

    double a = p4[0];
    double b = p7[1];

    uvCoords.resize(9);

    uvCoords[0][0] = -1.0;
    uvCoords[0][1] = -1.0;
    uvCoords[1][0] =  1.0;
    uvCoords[1][1] = -1.0;
    uvCoords[2][0] =  1.0;
    uvCoords[2][1] =  1.0;
    uvCoords[3][0] = -1.0;
    uvCoords[3][1] =  1.0;

    uvCoords[4][0] =  a;
    uvCoords[4][1] = -1.0;
    uvCoords[5][0] =  1.0;
    uvCoords[5][1] = b;
    uvCoords[6][0] =  a;
    uvCoords[6][1] =  1.0;
    uvCoords[7][0] =  -1.0;
    uvCoords[7][1] =  b;
    uvCoords[8][0] =   a;
    uvCoords[8][1] =  b;

    quads.resize(3);
    quads[0][0] = 1;
    quads[0][1] = 5;
    quads[0][2] = 8;
    quads[0][3] = 4;
    quads[1][0] = 2;
    quads[1][1] = 6;
    quads[1][2] = 8;
    quads[1][3] = 5;
    quads[2][0] = 3;
    quads[2][1] = 7;
    quads[2][2] = 8;
    quads[2][3] = 6;

    return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////

int JFEM2D :: splitSelfTangledCorner1( const JFacePtr &face, vector<Point2D> &uvCoords, vector<Array4I> &quads)
{
    uvCoords.clear();

    if( face->getSize(0) != 4 ) {
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
    zeroSearch( face, p0, p1, p4, found);
    if( !found ) return 2;

    p0[0] =   1.0;
    p0[1] =  -1.0;
    p1[0] =   1.0;
    p1[1] =   1.0;
    found = 0;
    Point2D p5;
    zeroSearch( face, p0, p1, p5, found);
    if( !found ) return 2;

    uvCoords.resize(9);

    double a = p4[0];
    double b = p5[1];
    uvCoords.resize(9);

    uvCoords[0][0] = -1.0;
    uvCoords[0][1] = -1.0;
    uvCoords[1][0] =  1.0;
    uvCoords[1][1] = -1.0;
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
    uvCoords[8][0] = a;
    uvCoords[8][1] =  b;

    quads.resize(3);
    quads[0][0] = 0;
    quads[0][1] = 4;
    quads[0][2] = 8;
    quads[0][3] = 7;

    quads[1][0] = 2;
    quads[1][1] = 6;
    quads[1][2] = 8;
    quads[1][3] = 5;

    quads[2][0] = 3;
    quads[2][1] = 7;
    quads[2][2] = 8;
    quads[2][3] = 6;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////

int JFEM2D :: splitSelfTangledCorner2( const JFacePtr &face, vector<Point2D> &uvCoords, vector<Array4I> &quads)
{
    uvCoords.clear();

    if( face->getSize(0) != 4 ) {
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
    zeroSearch( face, p0, p1, p5, found);
    if( !found ) return 2;

    p0[0] =  1.0;
    p0[1] =  1.0;
    p1[0] =  -1.0;
    p1[1] =   1.0;
    found = 0;
    Point2D p6;
    zeroSearch( face, p0, p1, p6, found);
    if( !found ) return 2;

    double a, b;
    a = p6[0];
    b = p5[1];

    uvCoords.resize(9);

    uvCoords[0][0] = -1.0;
    uvCoords[0][1] = -1.0;
    uvCoords[1][0] =  1.0;
    uvCoords[1][1] = -1.0;
    uvCoords[2][0] =  1.0;
    uvCoords[2][1] =  1.0;
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
    uvCoords[8][0] =  a;
    uvCoords[8][1] =  b;

    quads.resize(3);
    quads[0][0] = 0;
    quads[0][1] = 4;
    quads[0][2] = 8;
    quads[0][3] = 7;

    quads[1][0] = 1;
    quads[1][1] = 5;
    quads[1][2] = 8;
    quads[1][3] = 4;

    quads[2][0] = 3;
    quads[2][1] = 7;
    quads[2][2] = 8;
    quads[2][3] = 6;

    return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////

int JFEM2D :: splitSelfTangledCorner3( const JFacePtr &face, vector<Point2D> &uvCoords, vector<Array4I> &quads)
{
    uvCoords.clear();

    if( face->getSize(0) != 4 ) {
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
    zeroSearch( face, p0, p1, p6, found);
    if( !found ) return 2;

    p0[0] =  -1.0;
    p0[1] =   1.0;
    p1[0] =  -1.0;
    p1[1] =  -1.0;
    found = 0;
    Point2D p7;
    zeroSearch( face, p0, p1, p7, found);
    if( !found ) return 2;

    uvCoords.resize(9);

    double a = p6[0];
    double b = p7[1];

    uvCoords[0][0] = -1.0;
    uvCoords[0][1] = -1.0;
    uvCoords[1][0] =  1.0;
    uvCoords[1][1] = -1.0;
    uvCoords[2][0] =  1.0;
    uvCoords[2][1] =  1.0;
    uvCoords[3][0] = -1.0;
    uvCoords[3][1] =  1.0;

    uvCoords[4][0] =  a;
    uvCoords[4][1] = -1.0;
    uvCoords[5][0] = 1.0 ;
    uvCoords[5][1] = b;
    uvCoords[6][0] =  a;
    uvCoords[6][1] =  1.0;
    uvCoords[7][0] = -1.0;
    uvCoords[7][1] =  b;
    uvCoords[8][0] =  a;
    uvCoords[8][1] =  b;

    quads.resize(3);
    quads[0][0] = 0;
    quads[0][1] = 4;
    quads[0][2] = 8;
    quads[0][3] = 7;

    quads[1][0] = 1;
    quads[1][1] = 5;
    quads[1][2] = 8;
    quads[1][3] = 4;

    quads[2][0] = 2;
    quads[2][1] = 6;
    quads[2][2] = 8;
    quads[2][3] = 5;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////

int JFEM2D :: splitSelfTangledQuad( const JFacePtr &face, vector<Point2D> &uvCoords, vector<Array4I> &quads)
{
    Point2D uv;

    uv[0] = -1.0;
    uv[1] = -1.0;
    double J0 = getJacobian(face, uv);
    if( J0 < 0.0) return splitSelfTangledCorner0( face, uvCoords, quads);

    uv[0] = 1.0;
    uv[1] = -1.0;
    double J1 = getJacobian(face, uv);
    if( J1 < 0.0) return splitSelfTangledCorner1( face, uvCoords, quads);

    uv[0] = 1.0;
    uv[1] = 1.0;
    double J2 = getJacobian(face, uv);
    if( J2 < 0.0) return splitSelfTangledCorner2( face, uvCoords, quads);

    uv[0] =-1.0;
    uv[1] = 1.0;
    double J3 = getJacobian(face, uv);
    if( J3 < 0.0) return splitSelfTangledCorner3( face, uvCoords, quads);

    return 1;
}

void SelfTangledIntegrationRegions( const JFacePtr &face, JFaceSequence &faces)
{
}
#endif



