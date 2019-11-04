#include "FEM2D.hpp"
#include "Poisson2D.hpp"
#include "TangleElements2D.hpp"


////////////////////////////////////////////////////////////////////////////////
JPoisson2D::JPoisson2D()
{
    mesh   = nullptr;
    dofPerNode   = 1;
    setOrder(1);
}

////////////////////////////////////////////////////////////////////////////////

void JPoisson2D :: setDirichletValue( const JNodePtr &vtx, double u)
{
    int id = vtx->getID();
    dirichlet[id] = u;
}
////////////////////////////////////////////////////////////////////////////////

void JPoisson2D :: setDirichletValue( const JEdgePtr &edge, double u)
{
    edge->setAttribute("Dirichlet", u);
}

////////////////////////////////////////////////////////////////////////////////
void JPoisson2D :: setDirichletValue( const JEdgeSequence &vedges, double u)
{
    for( const JEdgePtr &e: vedges)  setDirichletValue(e, u);
}

////////////////////////////////////////////////////////////////////////////////
void JPoisson2D :: setNeumannValue( const JNodePtr &vtx, double u)
{
    int id = vtx->getID();
    neumann[id] = u;
}
////////////////////////////////////////////////////////////////////////////////

void JPoisson2D :: setNeumannValue( const JEdgePtr &edge, double u)
{
    edge->setAttribute("Neumann", u);
}

////////////////////////////////////////////////////////////////////////////////
void JPoisson2D :: setNeumannValue( const JEdgeSequence &vedges, double u)
{
    for( const JEdgePtr &e: vedges)  setNeumannValue(e, u);
}

////////////////////////////////////////////////////////////////////////////////

void JPoisson2D :: applyDirichletBC()
{
    int    nodeID;
    size_t numedges = mesh->getSize(1);

    JNodePtr vtx;

    double val = 0.0;
    for( size_t i = 0; i < numedges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        int err = edge->getAttribute("Dirichlet", val);
        if( !err) {
            nodeID = edge->getNodeAt(0)->getID();
            dirichlet[nodeID] = val;
            nodeID = edge->getNodeAt(1)->getID();
            dirichlet[nodeID] = val;
            if( shapeOrder == 2) {
                edge->getAttribute("Order2Node", vtx);
                nodeID  = vtx->getID();
                dirichlet[nodeID] = val;
            }
        }
    }

    if( dirichlet.empty() )
        cout << "Warning: There are no restricted nodes in the mesh " << endl;
}
////////////////////////////////////////////////////////////////////////////////

void JPoisson2D :: applyNeumannBC()
{
    if( verbose )
        cout << "Info: Assplying Neumann boundary condition " << endl;

    vector<double>  fb;
    vector<int>     globalPos;

    double val = 0.0;
    JEdgeSequence boundedges;
    mesh->getEntities("Neumann", boundedges);

    for( const JEdgePtr &edge : boundedges) {
        edge->getAttribute("Neumann", val);
        integrate(edge, val, fb);
        getAssemblyIndices(edge, globalPos);
        for( size_t i = 0; i < fb.size(); i++) {
            int gid  = globalPos[i];
            f[gid]  -= fb[i];
        }
    }

#ifdef DEBUG
    cout << " F-Vector: After Neumann BC" << endl;
    latexVector(f);
#endif

    if( verbose )
        cout << "Info: Boundary G Field applied condition applied" << endl;

}

/////////////////////////////////////////////////////////////////////////////////
void JPoisson2D :: applyBC()
{
    // Apply boundary force conditions ...
    applyNeumannBC();

    // Specify the fixed boundary condition....
    applyDirichletBC();
}
/////////////////////////////////////////////////////////////////////////////////

void JPoisson2D :: integrate( const JEdgePtr &edge, double value, vector<double> &fb)
{
    // Integrate "force" vector which is on the right hand side of the system...
    // Remember that Guass points are assigned from the left point to the right
    // point in the interval (-1, 1). But as per convention, the end nodes are
    // assembled first and then the nodes on the edge, therefore, we need to
    // make adjustment. For this reason, we first use "ftmp" and then assemble
    // it in proper order in the vector "fb"...
    int nnodes = numNodesPerElement[0];
    int ngauss = edgeGaussPoints.size();

    cout << nnodes << "  " << ngauss << endl;

    vector<double> ftmp(nnodes);
    for( int i = 0; i < nnodes; i++) ftmp[i] = 0.0;

    // For linear edges, the dJ can be calculated simply....
    double len = JEdgeGeometry::getLength( edge);
    double dJ = 0.5*len;

    // For each gauss point, evaluate the function f ( which is constant in
    // ou modelling).
    vector<double> N;
    for( int ig = 0; ig < ngauss; ig++) {
        double u   = edgeGaussPoints[ig][0];
        double w   = edgeGaussPoints[ig][1];
        double multby = w*value*dJ;
        femSpace->edgeElement.getShapeFunc(u, N);
        for( int ip = 0; ip < nnodes; ip++)
            ftmp[ip] +=  multby*N[ip];
    }

    // Rearrange now. First, two nodes of the edge and then higher order
    // nodes...
    fb.resize( nnodes );

    fb[0] = ftmp[0];
    fb[1] = ftmp[nnodes-1];
    for( int i = 1; i < nnodes-1; i++)
        fb[i+1] = ftmp[i];
}

/////////////////////////////////////////////////////////////////////////////////

void JPoisson2D :: getField( vector<double> &unode)
{
    unode.clear();
    if( mesh == nullptr || u.empty() ) return;

    int numnodes = mesh->getSize(0);
    unode.resize(numnodes);
    for( int i = 0; i < numnodes; i++)
        unode[i] = u[i];
}
////////////////////////////////////////////////////////////////////////////////////

void JPoisson2D :: BMatrix(const JFacePtr &face, const Point2D &uv, MatrixXd &B, double &dJ)
{
    MatrixXd gradN, xy;

    getShapeDeriv(face, uv, gradN);
    getCoords(face, xy);
    Matrix2d J  = gradN*xy;
    dJ =  J(0,0)*J(1,1) - J(0,1)*J(1,0);
    Matrix2d invJ = J.inverse();

    //------------------------------------------------------------------------------
    // Calculate :
    //   dNx_1   dNx_2  .....   dNx_p
    //   dNy_1   dNy_2  ....    dNy_p
    //------------------------------------------------------------------------------

    // Calculate the "B" Matrix at the Guass point on a face ...
    double j00 = invJ.coeffRef(0,0);
    double j01 = invJ.coeffRef(0,1);
    double j10 = invJ.coeffRef(1,0);
    double j11 = invJ.coeffRef(1,1);

    B.setZero();
    int nf = xy.rows();
    for( int i= 0; i < nf; i++) {
        double dNdU = gradN.coeffRef(0,i);   // Partial N over Partial xi
        double dNdV = gradN.coeffRef(1,i);   // Partial N over Partial eta
        double dNdX = j00*dNdU + j01*dNdV;    // partial N over Partial X
        double dNdY = j10*dNdU + j11*dNdV;    // partial N over Partial Y
        B.coeffRef(0,i) = dNdX;
        B.coeffRef(1,i) = dNdY;
    }
}

/////////////////////////////////////////////////////////////////////////////////////

void JPoisson2D :: BMatrix(const JFacePtr &face, const Point2D &uv, MatrixXd &B)
{
    double dJ;
    BMatrix(face, uv, B, dJ);
}

/////////////////////////////////////////////////////////////////////////////////////

void JPoisson2D :: assembleTangledK( const JFacePtr &face1, const JFacePtr &face2)
{
    assert( face1 != face2);

    int ndof = 0;
    if( face1->getSize(0) == 3)
        ndof = dofPerNode*numNodesPerElement[1];
    else
        ndof = dofPerNode*numNodesPerElement[2];

    MatrixXd Kp(ndof, ndof);
    Kp.setZero();

    int err = integrate(face1, face2, Kp);
    if( err ) return;

    vector<int> adof, bdof;
    getAssemblyIndices(face1, adof);
    getAssemblyIndices(face2, bdof);
    for( int i = 0; i < ndof; i++) {
        for( int j = 0; j < ndof; j++) {
            Ktangle.coeffRef(adof[i],bdof[j]) += Kp(i,j);
            Ktangle.coeffRef(bdof[i],adof[j]) += Kp(j,i);
        }
    }
}
////////////////////////////////////////////////////////////////////////////////////

void JPoisson2D :: assembleTangledK()
{
    Ktangle.resize(numDOF, numDOF);
    Ktangle.setZero();

#ifdef CSV
    /*
      tangledK = 0;
      numTangledPairs = 0;
      for( const std::pair<int,int> &boxPair: boxPairs) {
            const JFacePtr &face1 = mesh->getFaceAt(boxPair.first);
            const JFacePtr &face2 = mesh->getFaceAt(boxPair.second);
            assembleTangledK(face1, face2);
       }
    */
#endif

    if( verbose)
        cout << "Info: #of Tangled Pairs " << numTangledPairs << endl;
}
////////////////////////////////////////////////////////////////////////////////////


int JPoisson2D::solveLinearSystem()
{
    if( verbose ) {
        cout << "Info: Solving linear system ... " << endl;
        cout << "Info: #Nonzeros in Sparse Matrix " << K.nonZeros() << endl;
    }

    // For small size problems, we can use Direct solver ....
    vector<double> ureduced;
    int err =  LinearSystem::solve(K, freduced, ureduced);

    // Reassemble the vector "U" in full order (i.e. reduced + fixed value)..
    vector<double> ufull(numDOF);
    for( size_t i = 0; i < freeDof.size(); i++)
        ufull[freeDof[i]] = ureduced[i];

    for( const pair<int,double> &keyVal: dirichlet)
        ufull[keyVal.first] = keyVal.second;

    // Sperate "u" and "v" into indepedent vectors, This is what matters
    // in the final stage...
    size_t numNodes = mesh->getSize(0);
    numNodes +=  (numNodesPerElement[0]-2)*mesh->getSize(1);

    u.resize( numNodes);
    for( size_t i = 0; i < numNodes; i++) {
        u[i] = ufull[i];
    }

    if( verbose )
        cout << "Info: Linear system solved ... " << endl;

    return err;
}

/////////////////////////////////////////////////////////////////////////////////////
void JPoisson2D :: integrate_self_tangle(const JFacePtr &face, MatrixXd &Ke, vector<double> &fe)
{
    /*
       vector<Point2D>   uvCoords;
       vector<Array4I>   quads;

       splitSelfTangledQuad( face, uvCoords, quads);

       vector<Point3D>  gaussPoints(12);

       JNodeSequence nodes(4);
       for( int i = 0; i < 4; i++) nodes[i] = JVertex::newObject();
       JFacePtr qface = Quadrilateral::newObject(nodes);

       Point3D xyz;
       xyz[0] = 0.0;
       xyz[1] = 0.0;
       xyz[2] = 0.0;
       Point2D xy;
       Point2D uv;

       femSpace->getFaceGaussPoints(FEMSpace::QUAD, quadGaussPoints);

       int index = 0;
       for( int i = 0; i < 3; i++) {
            for( int j = 0; j < 4; j++) {
                 int id = quads[i][j];
                 xyz[0] = uvCoords[id][0];
                 xyz[1] = uvCoords[id][1];
                 qface->getNodeAt(j)->setXYZCoords(xyz);
            }
            for( int j = 0; j < 4; j++) {
                 uv[0] = quadGaussPoints[j][0];
                 uv[1] = quadGaussPoints[j][1];
                 getXYCoords(qface, uv, xy);
                 gaussPoints[index][0] = xy[0];
                 gaussPoints[index][1] = xy[1];
                 index++;
            }
        }

       vector<Point2D>   xyCoords(9);
       for( int i = 0; i < 9; i++)
            getXYCoords( face, uvCoords[i], xyCoords[i] );

       index = 0;
       for( int i = 0; i < 3; i++) {
            for( int j = 0; j < 4; j++) {
                 int id = quads[i][j];
                 xyz[0] = xyCoords[id][0];
                 xyz[1] = xyCoords[id][1];
                 qface->getNodeAt(j)->setXYZCoords(xyz);
            }
            for( int j = 0; j < 4; j++) {
                 uv[0] = quadGaussPoints[j][0];
                 uv[1] = quadGaussPoints[j][1];
                 double dJ  = getJacobian(qface, uv);
                 gaussPoints[index++][2] =  dJ;
            }
        }

        double sum = 0.0;
        for( int i = 0; i < 12; i++)
             sum += gaussPoints[i][2];

       cout << JFaceGeometry::getArea2D( face )  << " SUM = " << sum << endl;
        cout << "Done Upto hear" << endl;
        exit(0);
    */

}
/////////////////////////////////////////////////////////////////////////////////////


void JPoisson2D :: integrate(const JFacePtr &face, MatrixXd &Ke, vector<double> &fe)
{
    vector<Point3D> gaussPoints;

    int ndof, nnodes;
    if( face->getSize(0) == 3 ) {
        gaussPoints = triGaussPoints;
        nnodes      = numNodesPerElement[1];
    }

    if( face->getSize(0) == 4 ) {
        gaussPoints = quadGaussPoints;
        nnodes      = numNodesPerElement[2];
    }

    ndof = dofPerNode*nnodes;

    int ngauss = gaussPoints.size();

    //
    // Integrate over the face and return Element stiffness matrix "Ke" and "f" vector.
    // The size of matrix and f depends on the number of nodes.

    vector<double> dJ( ngauss);
    MatrixXd B;
    B.resize(2, ndof);

    Ke.resize( ndof, ndof);
    Ke.setZero();

    // Calculate the integral of (B'D*B*dA) on each face. The function
    // value is evaluated at the "Gauss Points". In our case, the
    // Inverse Jacobian is constant.

    Point2D uv;
    for( int ig = 0; ig < ngauss; ig++) {
        uv[0] = gaussPoints[ig][0];
        uv[1] = gaussPoints[ig][1];
        BMatrix( face, uv, B, dJ[ig]);
        if( absJacobian) dJ[ig] = fabs(dJ[ig]);
        double C = gaussPoints[ig][2]*dJ[ig];
        Ke = Ke + C*B.transpose()*B;
    }

    // Calculate the "f" vector, which depends on the body force ...
    fe.resize(ndof);
    for( int i = 0; i < ndof; i++)
        fe[i] = 0.0;

    // In classical method, we need to take +ve Jacobian. But when we have
    // tangled mesh, Jacobian will indicate sign of the element.
    /*
        if( !tangledK ) {
            for( int ig = 0; ig < ngauss; ig++)
                dJ[ig] = fabs(dJ[ig]);
        }
    */

    vector<double> N;
    Point2D xy;
    for( int ig = 0; ig < ngauss; ig++) {
        uv[0] = gaussPoints[ig][0];
        uv[1] = gaussPoints[ig][1];
        getShapeFunc(face, uv, N);
        getXYCoords( face, uv, xy);
        double val = fieldFunction->getScalar(xy);
        for( int ip = 0; ip < nnodes; ip++) {
            fe[ip] += gaussPoints[ig][2]*N[ip]*val*dJ[ig];
        }
    }

#ifdef DEBUG
    cout << "Calculating local stiffness matrix " << endl;
    latexMatrix(Ke);
    latexVector(fe);
#endif
}

/////////////////////////////////////////////////////////////////////////////////////

void JPoisson2D :: assembleK()
{
    if( mesh == nullptr) return;

    if( verbose)
        cout << "Assemble classical stiffness matrix starts ... " << endl;

    // Global matrix ...
    K.resize(numDOF, numDOF);
    K.setZero();

    MatrixXd   KElem;
    vector<int>   globalPos;
    vector<double> fElem;

    //
    // for each face, calculate the Ke and f and assemble them into global
    // matrix and vector respectively....
    //
    size_t numfaces = mesh->getSize(2);
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        integrate(face, KElem, fElem);
        getAssemblyIndices(face, globalPos);
        int ndof = globalPos.size();
        for( int j = 0; j < ndof; j++) {
            int irow = globalPos[j];
            for( int k = 0; k < ndof; k++) {
                int icol = globalPos[k];
                K.coeffRef(irow,icol) += KElem(j,k);
            }
            f[irow] += fElem[j];
        }
    }

#ifdef DEBUG
    cout << "Global stiffness matrix " << endl;
    latexMatrix(K);
    cout << "Force vector " << endl;
    latexVector(f);
#endif

    if( verbose)
        cout << "Stiffness matrix assembled, " << endl;
}

////////////////////////////////////////////////////////////////////////////////////
void JPoisson2D :: integrate( const JFacePtr &face1, const JFacePtr &face2,
                              const vector<Point2D> &triPnts, MatrixXd &Ke)
{
    /*
        assert( face1 != face2);
        assert( triPnts.size() == 3);

        double  dJ;
        Matrix2d invJ;
        // Find the area of subtriangle of the region...
        invJ(0,0)  = triPnts[2][1] - triPnts[0][1];
        invJ(0,1)  = triPnts[0][1] - triPnts[1][1];
        invJ(1,0)  = triPnts[0][0] - triPnts[2][0];
        invJ(1,1)  = triPnts[1][0] - triPnts[0][0];
        dJ   = invJ(0,0)*invJ(1,1) - invJ(0,1)*invJ(1,0);
        invJ = invJ/dJ;

        Ke.resize( dofPerFace, dofPerFace);
        Ke.setZero();

        if( fabs(dJ) < 1.0E-10) return;

        // Determine the orientation of the two faces ...
        double sign1 = JFaceGeometry::getOrientation2D( face1 );
        double sign2 = JFaceGeometry::getOrientation2D( face2 );
        assert(dJ > 0.0);

        // Determine the XY Coordinates of the gauss points in the triangle
        vector<Point2D> xy;
        tangleSpace->getFaceGaussXY(triPnts, xy);

        vector<double> u, v, weight;
        tangleSpace->getFaceGaussPoints( u, v, weight);

        Point2D uv;
        MatrixXd gradN, Bi, Bj;

        Bi.resize(3, dofPerFace);
        Bj.resize(3, dofPerFace);

        // Calculate the (u,v) of each gauss points in the two triangles
        // and calculate the Bi and Bj  Matrix..
        double coeff;
        int    ngauss = xy.size();

        for ( int ig = 0; ig < ngauss; ig++) {
            if( !isInside( triPnts, xy[ig]) ) {
                cout << "Fatal error : Gauss point not inside the triangle " << endl;
                exit(0);
            }
            getUVCoords( face1, xy[ig], uv );
            BMatrix(face1, uv, Bi);

            getUVCoords( face2, xy[ig], uv );
            BMatrix(face2, uv, Bj);

            coeff = sign1*sign2*weight[ig]*fabs(dJ);
            Ke += coeff*Bi.transpose()*Bj;
        }
    */
}

////////////////////////////////////////////////////////////////////////////////////

int JPoisson2D :: integrate( const JFacePtr &face1, const JFacePtr &face2, MatrixXd &Kp)
{
    /*
        assert( face1 != face2);

        Kp.setZero();

        // Check for intersection of the two faces ..
        vector<Point2D> polyPoints;
        getIntersection( face1, face2, polyPoints);

        int np = polyPoints.size();

        // If there was no intersection, no polygon will be formed...
        if( np < 3 ) return 1;

        // If the polygonal region is a triangle, no need to subdivide
        if( np == 3 ) {
            integrate(face1, face2, polyPoints, Kp);
            return 0;
        }

        // Since the intersection of two convex faces ( here triangles) is a also
        // a convex, we can simply triangulate the region by inserting a new vertex
        // at the center.

        Point2D center;
        center[0] = 0.0;
        center[1] = 0.0;
        for( int i = 0; i < np ; i++) {
            center[0] += polyPoints[i][0];
            center[1] += polyPoints[i][1];
        }
        center[0] /= (double)np;
        center[1] /= (double)np;

        assert( isInside(polyPoints, center));

        vector<Point2D> nodeCoords(3);

        // Integrate in the region..

        MatrixXd Ke(dofPerFace, dofPerFace);
        for( int i = 0; i < np; i++) {
            nodeCoords[0] = center;
            nodeCoords[1] = polyPoints[i];
            nodeCoords[2] = polyPoints[(i+1)%np];
            integrate(face1, face2, nodeCoords, Ke);
            Kp += Ke;
        }
    */
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JPoisson2D :: buildLinearSystem()
{
    if( mesh == nullptr ) return 1;

    K.setZero();
    Ktangle.setZero();

    f.resize(numDOF);

    for( int i = 0; i < numDOF; i++) f[i] = 0.0;

//  First calculate the Tangled K, so that we know whether we need to modify the
//  right hand side vector "f" in the classical formulation.
    assembleTangledK();


//  Classical assembly for tangle free mesh...
    assembleK();

    //
    // If the mesh is tangled, we will modify the K Matrix, Only the
    // internal edges are tangled and no edge must cross the boundary
    // edges ..
    //
    K = K + Ktangle;

    applyBC();


    // If there are large number of fixed value of "U" then we should reduce the
    // System for two reasons:
    //
    // (1) Run time performanec will increase with the reduced matrix size.
    // (2) Condition number will improve ...
    //
    freduced = f;
    reduceSystem(K, freduced, dirichlet, freeDof);

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
int JPoisson2D :: solve()
{
    return JFEM2D::solve();
}

///////////////////////////////////////////////////////////////////////////////
void JPoisson2D:: saveAs( const string &filename, const string &var)
{
    if( mesh == nullptr ) return;
    if( u.empty() ) return;

    size_t numnodes = mesh->getSize(0);

    double val;
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        vtx->setAttribute(var, u[i]);
    }

    JMeshVTKExporter mexp;
    mexp.addNodeAttribute(var);
    mexp.writeFile(mesh, filename);
}

/////////////////////////////////////////////////////////////////////////////////////
