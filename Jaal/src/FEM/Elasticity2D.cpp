#include "Elasticity2D.hpp"
#include "FEM2D.hpp"

int saveMatrix = 0;

#ifdef CSV
/////////////////////////////////////////////////////////////////////////////////////////
//
JElasticity2D::JElasticity2D()
{
    /*
        mesh   = nullptr;
        problem = PLAIN_STRESS;
        E    = 1.0E11;
        nu   = 0.33;
        F[0] = 0.0;
        F[1] = 0.0;
        setOrder(1);
        computeDMatrix();
        trishape.setOrder(1);
        bodyForce[0] = 0.0;
        bodyForce[1] = 0.0;
    */
}

///////////////////////////////////////////////////////////////////////////////

void JElasticity2D :: computeDMatrix()
{
    double C;

    if( problem == PLAIN_STRAIN ) {
        C =   E/((1+nu)*(1.0-2*nu));
        D(0,0) =   C*(1-nu);
        D(0,1) =   C*nu;
        D(0,2) =   0.0;

        D(1,0) =   C*nu;
        D(1,1) =   C*(1.0-nu);
        D(1,2) =   0.0;

        D(2,0) =   0.0;
        D(2,1) =   0.0;
        D(2,2) =   0.5*C*(1.0-2.0*nu);
    }

    if( problem == PLAIN_STRESS) {
        C =   E/(1-nu*nu);
        D(0,0) =   C;
        D(0,1) =   C*nu;
        D(0,2) =   0.0;

        D(1,0) =   C*nu;
        D(1,1) =   C;
        D(1,2) =   0.0;

        D(2,0) =   0.0;
        D(2,1) =   0.0;
        D(2,2) =   0.5*C*(1.0-nu);
    }

    if( debugStep == 1) {
        cout << "D Matrix " << endl;
        latexMatrix(D);
    }
}

/////////////////////////////////////////////////////////////////////////////////
void JElasticity2D :: initParams()
{
    /*
        nodesPerEdge       = femshape.getNumNodesPerEdge();
        nodesPerFace       = femshape.getNumNodesPerFace();

        dofPerEdge         = dofPerNode*nodesPerEdge;
        dofPerFace         = dofPerNode*nodesPerFace;
        numFaceGaussPoints = femshape.getNumFaceGaussPoints();
        numEdgeGaussPoints = femshape.getNumEdgeGaussPoints();
        femshape.getGaussPoints(exi, eweight);
        femshape.getGaussPoints(fxi, feta, fweight);

        u.clear();
        v.clear();
        nodeVonMises.clear();
        faceVonMises.clear();
        computeDMatrix();
    */
}

/////////////////////////////////////////////////////////////////////////////////

void JElasticity2D :: getDisplacements( vector<double> &unode, vector<double> &vnode)
{

    unode.clear();
    vnode.clear();
    if( mesh == nullptr || u.empty() || v.empty() ) return;

    int numnodes = mesh->getSize(0);
    unode.resize(numnodes);
    vnode.resize(numnodes);
    for( int i = 0; i < numnodes; i++) {
        unode[i] = u[i];
        vnode[i] = v[i];
    }
}

////////////////////////////////////////////////////////////////////////////////

void JElasticity2D :: setXForce( const JEdgeSequence &boundaryEdges, double force)
{
    // Specify the force on "geometric edge". The mesh edges on a geometric edge
    // gets uniform force.

    size_t nSize = boundaryEdges.size();

    if( nSize < 1) {
        cout << "Info: X-force boundary edges are empty " << endl;
        return;
    }

    double len = 0.0;
    for( size_t i = 0; i < nSize; i++) {
        const JEdgePtr &edge = boundaryEdges[i];
        len += JEdgeGeometry::getLength( edge );
    }

    for( size_t i = 0; i < nSize; i++) {
        const JEdgePtr &edge = boundaryEdges[i];
        edge->setAttribute("Xforce", force/len);
    }

}

/////////////////////////////////////////////////////////////////////////////////////

void JElasticity2D :: setYForce( const JEdgeSequence &boundaryEdges, double force)
{
    // Specify the force on "geometric edge". The mesh edges on a geometric edge
    // gets uniform force.

    size_t nSize = boundaryEdges.size();
    if( nSize < 1) {
        cout << "Info: X-force boundary edges are empty " << endl;
        return;
    }

    double len = 0.0;
    for( size_t i = 0; i < nSize; i++) {
        const JEdgePtr &edge = boundaryEdges[i];
        len += JEdgeGeometry::getLength( edge );
    }

    for( size_t i = 0; i < nSize; i++) {
        const JEdgePtr &edge = boundaryEdges[i];
        edge->setAttribute("Yforce", force/len);
    }
}

/////////////////////////////////////////////////////////////////////////////////////

void JElasticity2D :: fixX( const JNodePtr &vtx, double u)
{
    int id = vtx->getID();
    fixMap[2*id] = u;
}
/////////////////////////////////////////////////////////////////////////////////////
void JElasticity2D :: fixY( const JNodePtr &vtx, double v)
{
    int id = vtx->getID();
    fixMap[2*id+1] = v;
}
/////////////////////////////////////////////////////////////////////////////////////

void JElasticity2D :: fixX( const JEdgeSequence &boundaryEdges, double u)
{
    // Fix the "U" component of a edge.
    size_t nSize = boundaryEdges.size();
    for( size_t i = 0; i < nSize; i++) {
        const JEdgePtr &edge = boundaryEdges[i];
        edge->setAttribute("Xfixed", u);
    }
}
/////////////////////////////////////////////////////////////////////////////////////
void JElasticity2D :: fixY( const JEdgeSequence &boundaryEdges, double v)
{
    // Fix the "V" component of a edge.
    size_t nSize = boundaryEdges.size();
    for( size_t i = 0; i < nSize; i++) {
        const JEdgePtr &edge = boundaryEdges[i];
        edge->setAttribute("Yfixed", v);
    }
}
/////////////////////////////////////////////////////////////////////////////////////

void JElasticity2D :: fixEdges( const JEdgeSequence &boundaryEdges, double uv)
{
    // Specify which edges are completely restricted in both the directions...
    size_t nSize = boundaryEdges.size();
    for( size_t i = 0; i < nSize; i++) {
        const JEdgePtr &edge = boundaryEdges[i];
        edge->setAttribute("Xfixed", uv);
        edge->setAttribute("Yfixed", uv);
    }
}

/////////////////////////////////////////////////////////////////////////////////////

void JElasticity2D :: BMatrix(const JFacePtr &face, const Point2D &uv, MatrixXd &B, double &dJ)
{
    /*
        int nsize = face->getSize(0);

        MatrixXd gradN;
        femshape.getShapeDeriv( uv[0], uv[1], gradN);

        MatrixXd xy(nsize,2);
        for( int i = 0; i < nsize; i++) {
            const Point3D &xyz = face->getNodeAt(i)->getXYZCoords();
            xy(i,0) = xyz[0];
            xy(i,1) = xyz[1];
        }

        Matrix2d J  = gradN.block(0,0,2,nsize)*xy;
        dJ =  J(0,0)*J(1,1) - J(0,1)*J(1,0);
        Matrix2d invJ = J.inverse();

        //------------------------------------------------------------------------------
        // Calculate :
        //   dNx_1   0     dNx_2   0         .....   dNx_p  0
        //   0     dNy_1   0       dNy_2     ....    0      dNy_p
        //   dNy_1 dNx_1   dNy_2   dNx_2     .....   dNy_p  dNx_p
        //------------------------------------------------------------------------------

        assert( gradN.rows() == 2);
        assert( gradN.cols() == nodesPerFace);

        // Calculate the "B" Matrix at the Guass point on a face ...
        double j00 = invJ.coeff(0,0);
        double j01 = invJ.coeff(0,1);
        double j10 = invJ.coeff(1,0);
        double j11 = invJ.coeff(1,1);

        B.setZero();
        for( int i= 0; i < nodesPerFace; i++) {
            double dNdU = gradN.coeff(0,i);   // Partial N over Partial xi
            double dNdV = gradN.coeff(1,i);   // Partial N over Partial eta
            double dNdX = j00*dNdU + j01*dNdV;    // partial N over Partial X
            double dNdY = j10*dNdU + j11*dNdV;    // partial N over Partial Y
            B.coeffRef(0,2*i)   = dNdX;
            B.coeffRef(1,2*i+1) = dNdY;
            B.coeffRef(2,2*i)   = dNdY;
            B.coeffRef(2,2*i+1) = dNdX;
        }
    */
}

/////////////////////////////////////////////////////////////////////////////////////

void JElasticity2D :: BMatrix(const JFacePtr &face, const Point2D &uv, MatrixXd &B)
{
    double dJ;
    BMatrix(face, uv, B, dJ);
}

/////////////////////////////////////////////////////////////////////////////////////

void JElasticity2D :: convexIntegrate(const JFacePtr &face, MatrixXd &Ke, vector<double> &fe)
{
    int err = 1;
    if( face->getSize(0) == 3 ) {
        if( numFaceGaussPoints == 1 || numFaceGaussPoints == 3  || numFaceGaussPoints == 7 ) err = 0;
    }

    if( face->getSize(0) == 4 ) {
        if( numFaceGaussPoints == 4 ) err = 0;
    }

    if( err ) {
        cout << " Warning: Invalid number of Gauss points for the element " <<endl;
        abort();
    }

    Ke.setZero();

    //
    // Integrate over the face and return Element stiffness matrix "Ke" and "f" vector.
    // The size of matrix and f depends on the number of nodes.

    vector<double> dJ( numFaceGaussPoints);

    MatrixXd B;
    B.resize(3, dofPerFace);

    // Calculate the integral of (B'D*B*dA) on each face. The function
    // value is evaluated at the "Gauss Points". In our case, the
    // Invere Jacobian is constant.

    Point2D uv;
    double  C;
    for( int ig = 0; ig < numFaceGaussPoints; ig++) {
        uv[0] = fxi[ig];
        uv[1] = feta[ig];
        BMatrix( face, uv, B, dJ[ig]);
        if( absJacobian )
            C = fweight[ig]*fabs(dJ[ig]);    // Commericial codes use fabs
        else
            C = fweight[ig]*dJ[ig];          // New observation: Solves many problems.
        Ke = Ke + C*B.transpose()*D*B;
    }

    // Calculate the "f" vector, which depends on the body force ...
    fe.resize(dofPerFace);
    for( int i = 0; i < dofPerFace; i++)
        fe[i] = 0.0;

    // In classical method, we need to take +ve Jacobian. But when we have
    // tangled mesh, Jacobian will indicate sign of the element.
    /*
        if( !tangledK ) {
            for( int ig = 0; ig < numFaceGaussPoints; ig++)
                dJ[ig] = fabs(dJ[ig]);
    //            dJ[ig] = dJ[ig];
        }
    */

    double bx = bodyForce[0];
    double by = bodyForce[1];
    for( int ip = 0; ip < nodesPerFace; ip++) {
        for( int ig = 0; ig < numFaceGaussPoints; ig++) {
            fe[2*ip]   += fweight[ig]*NFace[ig][ip]*bx*dJ[ig];
            fe[2*ip+1] += fweight[ig]*NFace[ig][ip]*by*dJ[ig];
        }
    }

    if( debugStep == 2 ) {
        cout << "Calculating local stiffness matrix " << endl;
        latexMatrix(Ke);
        latexVector(fe);
    }
}

//////////////////////////////////////////////////////////////////////////////////////////

void JElasticity2D :: concaveQuadTriIntegrate(const JFacePtr &quadface, int triangle, MatrixXd &Ke, vector<double> &fe)
{
    //
    // Integrate over the face and return Element stiffness matrix "Ke" and "f" vector.
    // The size of matrix and f depends on the number of nodes.

    /*
        trishape.setNumFaceGaussPoints(3);

        vector<double> ut, vt, weight;
        trishape.getGaussPoints(ut, vt, weight);

        int nGauss = ut.size();
        assert( nGauss == 3 );

        MatrixXd xyCoords(3,2), qParam(3,2);

        Point2D  xy;

        if( triangle == 0) {
            xy = quadface->getNodeAt(0)->getXYCoords();
            xyCoords(0,0) = xy[0];
            xyCoords(0,1) = xy[1];
            qParam(0,0)   = -1.0;
            qParam(0,1)   = -1.0;

            xy = quadface->getNodeAt(1)->getXYCoords();
            xyCoords(1,0) = xy[0];
            xyCoords(1,1) = xy[1];
            qParam(1,0)   =  1.0;
            qParam(1,1)   = -1.0;

            xy = quadface->getNodeAt(2)->getXYCoords();
            xyCoords(2,0) = xy[0];
            xyCoords(2,1) = xy[1];
            qParam(2,0)   =  1.0;
            qParam(2,1)   =  1.0;
        } else {
            xy = quadface->getNodeAt(0)->getXYCoords();
            xyCoords(0,0) = xy[0];
            xyCoords(0,1) = xy[1];
            qParam(0,0)   =  -1.0;
            qParam(0,1)   =  -1.0;

            xy = quadface->getNodeAt(2)->getXYCoords();
            xyCoords(1,0) = xy[0];
            xyCoords(1,1) = xy[1];
            qParam(1,0)   =  1.0;
            qParam(1,1)   =  1.0;

            xy = quadface->getNodeAt(3)->getXYCoords();
            xyCoords(2,0) = xy[0];
            xyCoords(2,1) = xy[1];
            qParam(2,0)   = -1.0;
            qParam(2,1)   =  1.0;
        }

        MatrixXd B;
        B.resize(3, dofPerFace);

        Matrix2d invJ;
        // Find the area of subtriangle of the region...
        invJ(0,0)  = xyCoords(2,1) - xyCoords(0,1);
        invJ(0,1)  = xyCoords(0,1) - xyCoords(1,1);
        invJ(1,0)  = xyCoords(0,0) - xyCoords(2,0);
        invJ(1,1)  = xyCoords(1,0) - xyCoords(0,0);
        double dJ  = invJ(0,0)*invJ(1,1) - invJ(0,1)*invJ(1,0);

        assert( fabs(dJ) > 0.0);

        for( int i = 0; i < 4; i++) {
            xy = quadface->getNodeAt(i)->getXYCoords();
        }

        // Calculate the integral of (B'D*B*dA) on each face. The function
        // value is evaluated at the "Gauss Points". In our case, the
        // Inverse Jacobian is constant.

        Point2D uvtri, uvquad;
        MatrixXd gradNT, gradNQ;
        Matrix2d  J1, J2, J1J2;
        B.setZero();
    */

//   assert(femshape.getElementType() == 4);

    /*
        double  C;
        for( int ig = 0; ig < nGauss; ig++) {
            uvtri[0] = ut[ig];
            uvtri[1] = vt[ig];
            trishape.getShapeDeriv( uvtri[0], uvtri[1], gradNT);
            assert( gradNT.rows() == 2);
            assert( gradNT.cols() == 3);
            J1 = gradNT*xyCoords;
            J2 = gradNT*qParam;
            assert( J1.rows() == 2);
            assert( J1.cols() == 2);
            mapUV( triangle, uvtri, uvquad);
            femshape.getShapeDeriv( uvquad[0], uvquad[1], gradNQ);
            assert( gradNQ.rows() == 2);
            assert( gradNQ.cols() == 4);
            J1J2 = J1.inverse()*J2;
            assert( J1J2.rows() == 2);
            assert( J1J2.cols() == 2);
            for( int i = 0; i < 4; i++) {
                double dNQdx = J1J2(0,0)*gradNQ(0,i) + J1J2(0,1)*gradNQ(1,i);
                double dNQdy = J1J2(1,0)*gradNQ(0,i) + J1J2(1,1)*gradNQ(1,i);
                B(0,2*i+0) = dNQdx;
                B(1,2*i+1) = dNQdy;
                B(2,2*i+0) = dNQdy;
                B(2,2*i+1) = dNQdx;
            }

            if( absJacobian )
                C = fweight[ig]*fabs(dJ);    // Commericial codes use fabs
            else
                C = fweight[ig]*dJ;          // New observation: Solves many problems.
            Ke = Ke + C*B.transpose()*D*B;
        }
    */


    /*  We need to update this when there are body forces ...
        if( !tangledK )
             dJ = fabs(dJ);
        }
        double bx = bodyForce[0];
        double by = bodyForce[1];
        for( int ip = 0; ip < nodesPerFace; ip++) {
            for( int ig = 0; ig < nGaussPnts; ig++) {
                fe[2*ip]   += fweight[ig]*NFace[ig][ip]*bx*dJ;
                fe[2*ip+1] += fweight[ig]*NFace[ig][ip]*by*dJ;
            }
        }
    */
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void JElasticity2D :: meanValueTriIntegrate(const JFacePtr &quadface, int triangle, MatrixXd &Ke, vector<double> &fe)
{
    /*
        int nGauss = 3;
        trishape.setNumFaceGaussPoints(nGauss);

        vector<double> xi, eta, fweight;
        trishape.getGaussPoints(xi, eta, fweight);

        vector<Point2D>  triCoords(3);
        Point2D  xy;
        if( triangle == 0) {
            xy = quadface->getNodeAt(0)->getXYCoords();
            triCoords[0][0] = xy[0];
            triCoords[0][1] = xy[1];

            xy = quadface->getNodeAt(1)->getXYCoords();
            triCoords[1][0] = xy[0];
            triCoords[1][1] = xy[1];

            xy = quadface->getNodeAt(2)->getXYCoords();
            triCoords[2][0] = xy[0];
            triCoords[2][1] = xy[1];
        } else {
            xy = quadface->getNodeAt(0)->getXYCoords();
            triCoords[0][0]  = xy[0];
            triCoords[0][1]  = xy[1];

            xy = quadface->getNodeAt(2)->getXYCoords();
            triCoords[1][0] = xy[0];
            triCoords[1][1] = xy[1];

            xy = quadface->getNodeAt(3)->getXYCoords();
            triCoords[2][0]  = xy[0];
            triCoords[2][1]  = xy[1];
        }

        Matrix2d invJ;
        // Find the area of subtriangle of the region...
        invJ(0,0)  = triCoords[2][1] - triCoords[0][1];
        invJ(0,1)  = triCoords[0][1] - triCoords[1][1];
        invJ(1,0)  = triCoords[0][0] - triCoords[2][0];
        invJ(1,1)  = triCoords[1][0] - triCoords[0][0];
        double dJ  = (invJ(0,0)*invJ(1,1) - invJ(0,1)*invJ(1,0));

        if( fabs(dJ) < 0.0) return;

        // Calculate the integral of (B'D*B*dA) on each face. The function
        // value is evaluated at the "Gauss Points". In our case, the
        // Inverse Jacobian is constant.

        vector<Point2D> samplePoints;
        trishape.getFaceGaussXY( triCoords, samplePoints);

        int nSize = quadface->getSize(0);
        assert( nSize == 4 );
        vector<Point2D>  polyCoords(nSize);
        for( int i = 0; i < nSize; i++)
            polyCoords[i] = quadface->getNodeAt(i)->getXYCoords();

        JMeanValueCoordinates mvc;

        MatrixXd B;
        B.resize(3, 2*nSize);
        B.setZero();

        vector<Vec2D> gradNT;
        double  C;
        for( int ig = 0; ig < nGauss; ig++) {
            int err = mvc.getShapeGradients( polyCoords, samplePoints[ig], gradNT);
            if( !err ) {
                for( int i = 0; i < nSize; i++) {
                    double dNQdx = gradNT[i][0];
                    double dNQdy = gradNT[i][1];
                    B(0,2*i+0) = dNQdx;
                    B(1,2*i+1) = dNQdy;
                    B(2,2*i+0) = dNQdy;
                    B(2,2*i+1) = dNQdx;
                }
                C = fweight[ig]*fabs(dJ);
                Ke = Ke + C*B.transpose()*D*B;
            }
        }
    */

    /*  We need to update this when there are body forces ...
        if( !tangledK )
             dJ = fabs(dJ);
        }
        double bx = bodyForce[0];
        double by = bodyForce[1];
        for( int ip = 0; ip < nodesPerFace; ip++) {
            for( int ig = 0; ig < nGaussPnts; ig++) {
                fe[2*ip]   += fweight[ig]*NFace[ig][ip]*bx*dJ;
                fe[2*ip+1] += fweight[ig]*NFace[ig][ip]*by*dJ;
            }
        }
    */
}


//////////////////////////////////////////////////////////////////////////////////////////

void JElasticity2D :: meanValueIntegrate(const JFacePtr &face, MatrixXd &Ke, vector<double> &fe)
{
    /*
        Ke.setZero();
        fe.resize(dofPerFace);
        for( int i = 0; i < dofPerFace; i++)
            fe[i] = 0.0;

        if( face->getSize(0) == 3 ) {
            meanValueTriIntegrate(face, 0, Ke, fe);
            return;
        }

        assert( face->getSize(0) == 4 );

        // Split the quad into two triangles: The concave corner becomes the divider...
        int corner = getConcaveCorner(face);
        if( corner >= 0) {
            const JNodePtr &p = face->getNodeAt(corner+2);
            face->setStartNode(p);
        }

        meanValueTriIntegrate(face, 0, Ke, fe);
        meanValueTriIntegrate(face, 1, Ke, fe);

        if( debugStep == 2 ) {
            cout << "Calculating local stiffness matrix " << endl;
            latexMatrix(Ke);
            latexVector(fe);
        }
    */
}

//////////////////////////////////////////////////////////////////////////////////////////

void JElasticity2D :: integrate(const JFacePtr &face, MatrixXd &Ke, vector<double> &fe)
{
    if( shape_family == 0) {
        convexIntegrate(face, Ke, fe);
        return;
    }

    if(  shape_family == 2 ) {
        meanValueIntegrate(face, Ke, fe);
        return;
    }

    if( isConvex(face) ) {
        convexIntegrate(face, Ke, fe);
    } else {
        meanValueIntegrate(face, Ke, fe);
    }
}
/////////////////////////////////////////////////////////////////////////////////////

void JElasticity2D :: assembleK()
{
    /*
        if( mesh == nullptr) return;

        if( verbose)
            cout << "Assemble classical stiffness matrix starts ... " << endl;

        stopWatch.reset();
        stopWatch.start();

        // Calculate the stiffness matrix of each face and assemble into a global matrix..
        int numFaces = mesh->getSize(2);

        // Since all the elements are of same shape the Guass point location and weights
        // are constant, so calcualte them once..
        femshape.getGaussPoints(fxi, feta, fweight);

        // Calculate the shape function (N) and its derivative in UV space at each Gauss
        // point...
        NFace.resize(numFaceGaussPnts);
        gradUV.resize(numFaceGaussPnts);
        for( int i = 0; i < numFaceGaussPnts; i++) {
            femshape.getShapeFunc(  fxi[i], feta[i], NFace[i]  );
            femshape.getShapeDeriv( fxi[i], feta[i], gradUV[i] );
        }

        // Global matrix ...
        K.resize(numDOF, numDOF);

        // Element matrix ( local );
        MatrixXd KElem(dofPerFace, dofPerFace);
        vector<int> globalPos;
        vector<double> fElem(dofPerFace);

        //
        // for each face, calculate the Ke and f and assemble them into global
        // matrix and vector respectively....
        //

        for( int i = 0; i < numFaces; i++) {
            const JFacePtr &face = mesh->getFaceAt(i);
            integrate(face, KElem, fElem);
            getAssemblyIndices(face, globalPos);
            for( int j = 0; j < dofPerFace; j++) {
                int irow = globalPos[j];
                for( int k = 0; k < dofPerFace; k++) {
                    int icol = globalPos[k];
                    K.coeffRef(irow,icol) += KElem(j,k);
                }
                f[irow] += fElem[j];
            }
        }

        if( debugStep == 3 ) {
            cout << "Global stiffness matrix " << endl;
            latexMatrix(K);
            cout << "Force vector " << endl;
            latexVector(f);
        }

        stopWatch.stop();
        if( verbose) {
            cout << "Stiffness matrix assembled, " << endl;
            cout << "Execution time in assembling classical K matrix(sec) "
                 << stopWatch.getSeconds() << endl;
        }
    */

}

////////////////////////////////////////////////////////////////////////////////////

void JElasticity2D :: integrate( const JEdgePtr &edge, double force, vector<double> &fb)
{
    // Integrate "force" vector which is on the right hand side of the system...
    // Remember that Guass points are assigned from the left point to the right
    // point in the interval (-1, 1). But as per convention, the end nodes are
    // assembled first and then the nodes on the edge, therefore, we need to
    // make adjustment. For this reason, we first use "ftmp" and then assemble
    // it in proper order in the vector "fb"...

    /*
        vector<double> ftmp(nodesPerEdge);
        for( int i = 0; i < nodesPerEdge; i++) ftmp[i] = 0.0;

        // For linear edges, the dJ can be calculated simply....
        double len = JEdgeGeometry::getLength( edge );
        double dJ = 0.5*len;

        // For each gauss point, evaluate the function f ( which is constant in
        // ou modelling).
        for( int ig = 0; ig < numEdgeGaussPnts; ig++) {
            for( int ip = 0; ip < nodesPerEdge; ip++)
                ftmp[ip] +=  eweight[ig]*NEdge[ig][ip]*force*dJ;
        }

        // Rearrange now. First, two nodes of the edge and then higher order
        // nodes...
        fb.resize( nodesPerEdge);

        fb[0] = ftmp[0];
        fb[1] = ftmp[nodesPerEdge-1];
        for( int i = 1; i < nodesPerEdge-1; i++)
            fb[i+1] = ftmp[i];
    */
}

////////////////////////////////////////////////////////////////////////////////////

void JElasticity2D :: applyNeumannBC()
{
    /*
        if( verbose )
            cout << "Info: Assplying force boundary condition, if specified " << endl;

        femshape.getGaussPoints(exi, eweight);

        NEdge.resize(numEdgeGaussPnts);
        gradU.resize(numEdgeGaussPnts);
        for( int i = 0; i < numEdgeGaussPnts; i++)
            femshape.getShapeFunc(  exi[i], NEdge[i]  );

        vector<double>  fb;
        vector<int>     globalPos;

        JEdgeSequence boundedges;

        int nCount = 0;
        double val = 0.0;
        mesh->getEntities("Xforce", boundedges);
        nCount += boundedges.size();
        for( const JEdgePtr &edge : boundedges) {
            int err = edge->getAttribute("Xforce", val);
            if( !err) {
                integrate(edge, val, fb);
                getAssemblyIndices(edge, globalPos);
                for( size_t i = 0; i < fb.size(); i++) {
                    int gid  = globalPos[2*i];
                    f[gid]  += fb[i];
                }
            }
        }

        mesh->getEntities("Yforce", boundedges);
        nCount += boundedges.size();
        for( const JEdgePtr &edge : boundedges) {
            int err = edge->getAttribute("Yforce", val);
            if( !err) {
                integrate( edge, val, fb);
                getAssemblyIndices(edge, globalPos);
                for( size_t i = 0; i < fb.size(); i++) {
                    int gid = globalPos[2*i+1];
                    f[gid] += fb[i];
                }
            }
        }

        if( debugStep == 4) {
            cout << " F-Vector: After Neumann BC" << endl;
            latexVector(f);
        }

        if( verbose && nCount )
            cout << "Info: Force boundary condition applied" << endl;
    */
}
////////////////////////////////////////////////////////////////////////////////////

void JElasticity2D :: applyDirichletBC()
{
    int    nodeID;
    size_t numedges = mesh->getSize(1);

    for( size_t i = 0; i < numedges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        if( edge->hasAttribute("Xfixed") ) {
            nodeID = edge->getNodeAt(0)->getID();
            fixMap[2*nodeID] = 0.0;
            nodeID = edge->getNodeAt(1)->getID();
            fixMap[2*nodeID] = 0.0;
            if( shapeOrder == FEMSpace::QUADRATIC) {
                int err = edge->getAttribute("Order2Node", nodeID);
                if( !err) fixMap[2*nodeID] = 0.0;
            }
        }
        if( edge->hasAttribute("Yfixed") ) {
            nodeID = edge->getNodeAt(0)->getID();
            fixMap[2*nodeID+1] = 0.0;
            nodeID = edge->getNodeAt(1)->getID();
            fixMap[2*nodeID+1] = 0.0;
            if( shapeOrder == FEMSpace::QUADRATIC) {
                int err = edge->getAttribute("Order2Node", nodeID);
                if( !err) fixMap[2*nodeID+1] = 0.0;
            }
        }
    }

    if( fixMap.empty() )
        cout << "Warning: There are no restricted nodes in the mesh " << endl;

    if( verbose)
        cout << "Info: Dirichlet boundary conditions applied  " << endl;
}

////////////////////////////////////////////////////////////////////////////////////

void JElasticity2D :: applyBC()
{
    if( verbose) {
        stopWatch.reset();
        stopWatch.start();
    }

    // Apply boundary force conditions ...
    applyNeumannBC();

    // Specify the fixed boundary condition....
    applyDirichletBC();

    if( verbose) {
        stopWatch.stop();
        cout << "Execution time in applying boundary conditions (sec) "
             << stopWatch.getSeconds() << endl;
    }
}

////////////////////////////////////////////////////////////////////////////////////
void JElasticity2D :: integrate( const JFacePtr &face1, const JFacePtr &face2,
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
        trishape.getFaceGaussXY(triPnts, xy);

        vector<double> u, v, weight;
        trishape.getGaussPoints( u, v, weight);

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
            Ke += coeff*Bi.transpose()*D*Bj;
        }
    */
}

////////////////////////////////////////////////////////////////////////////////////

int JElasticity2D :: integrate( const JFacePtr &face1, const JFacePtr &face2, MatrixXd &Kp)
{
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
        center[0] += 0.5*(polyPoints[(i+1)%np][0] + polyPoints[(i+2)%np][0]);
        center[1] += 0.5*(polyPoints[(i+1)%np][1] + polyPoints[(i+2)%np][1]);
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

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

void JElasticity2D :: assembleTangledK( const JFacePtr &face1, const JFacePtr &face2)
{
    /*
        assert( face1 != face2);

        MatrixXd Kp(dofPerFace, dofPerFace);

        Kp.setZero();

        assert(Kp.rows() == dofPerFace);
        assert(Kp.cols() == dofPerFace);

        int err = integrate(face1, face2, Kp);

        if( !err) {
            vector<int> adof, bdof;
            getAssemblyIndices(face1, adof);
            getAssemblyIndices(face2, bdof);
            for( int i = 0; i < dofPerFace; i++) {
                for( int j = 0; j < dofPerFace; j++) {
                    Ktangle.coeffRef(adof[i],bdof[j]) += Kp(i,j);
                    Ktangle.coeffRef(bdof[i],adof[j]) += Kp(j,i);
                }
            }
            tangledK = 1;
            numTangledPairs++;
        }
    */
}
////////////////////////////////////////////////////////////////////////////////////

void JElasticity2D :: assembleTangledK()
{

#ifdef CSV
    if( verbose) {
        stopWatch.reset();
        stopWatch.start();
    }

    tangledK = 0;
    numTangledPairs = 0;

    int numDOF =  dofPerNode*mesh->getSize(0);

    Ktangle.resize(numDOF, numDOF);
    Ktangle.setZero();
    cout << "Debug exit" << endl;
    exit(0);

    /*
        for( const std::pair<int,int> &boxPair: boxPairs) {
            const JFacePtr &face1 = mesh->getFaceAt(boxPair.first);
            const JFacePtr &face2 = mesh->getFaceAt(boxPair.second);
            assembleTangledK(face1, face2);
        }
    */

    if( verbose) {
        cout << "Info: #of Tangled Pairs " << numTangledPairs << endl;
        stopWatch.stop();
        cout << "Execution time in assembling tangled K matrix " << stopWatch.getSeconds() << endl;
    }
#endif
}
////////////////////////////////////////////////////////////////////////////////////
int JElasticity2D :: buildLinearSystem()
{
    if( mesh == nullptr ) return 1;

    initParams();

    int numNodes = mesh->getSize(0);
    int numEdges = mesh->getSize(1);
    assert( numNodes );
    assert( numEdges );

    numDOF = dofPerNode*(numNodes + numEdges*(nodesPerEdge-2));

    f.resize(numDOF);

    for( int i = 0; i < numDOF; i++) f[i] = 0.0;

    if( shapeOrder == FEMSpace::QUADRATIC) genQuadraticNodes();

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
    reduceSystem(K, freduced, fixMap, freeDof);
    return 0;
}

int JElasticity2D::solveLinearSystem()
{

    if( verbose ) {
        cout << "Info: Solving linear system ... " << endl;
        cout << "Info: #Nonzeros in Sparse Matrix " << K.nonZeros() << endl;
        stopWatch.reset();
        stopWatch.start();
    }

    static int  nCount = 0;

    if( saveMatrix ) {
        ostringstream oss;
        oss << "K";
        if( nCount < 10) oss << "0";
        oss << nCount++;
        oss << ".dat";
        ofstream ofile(oss.str().c_str(), ios::out);
        ofile << "# name: K" << endl;
        ofile << "# type: sparse matrix" << endl;
        ofile << "# nnz: " << K.nonZeros() << endl;
        ofile << "# rows: " << K.rows() << endl;
        ofile << "# columns: " << K.cols() << endl;
        for( int j = 0; j < K.cols(); j++) {
            for( int i = 0; i < K.rows(); i++) {
                double val = K.coeff(i,j);
                if( val ) ofile << i+1 << " " << j+1 << " " << val << endl;
            }
        }
    }
    // For small size problems, we can use Direct solver ....
    vector<double> ureduced;
    int err =  LinearSystem::solve( K, freduced, ureduced);

    // Reassemble the vector "U" in full order (i.e. reduced + fixed value)..
    vector<double> ufull(numDOF);
    for( size_t i = 0; i < freeDof.size(); i++)
        ufull[freeDof[i]] = ureduced[i];

    for( const pair<int,double> &keyVal: fixMap)
        ufull[keyVal.first] = keyVal.second;

    // Sperate "u" and "v" into indepedent vectors, This is what matters
    // in the final stage...
    int numNodes = mesh->getSize(0) + (nodesPerEdge-2)*mesh->getSize(1);
    u.resize( numNodes);
    v.resize( numNodes);
    for( int i = 0; i < numNodes; i++) {
        u[i] = ufull[2*i];
        v[i] = ufull[2*i+1];
    }

    if( debugStep == 7) {
        cout << "Displacement vector " << endl;
        latexVector(ufull);
    }

    if( verbose ) {
        cout << "Info: Linear system solved ... " << endl;
        stopWatch.stop();
        cout << "Execution time in solving linear equations " << stopWatch.getSeconds() << endl;
    }

    return err;
}

/////////////////////////////////////////////////////////////////////////////////////

void JElasticity2D:: computeStress( const JFacePtr &face)
{
    /*
        // Once the displacement are know, we can calculate the strain and stress
        // in each face. For Von-Mises stress is the maximum stress at any gauss point..
        int fid = face->getID();

        VectorXd un(2*nodesPerFace);
        vector<int> nodes;
        getNodes(face, nodes);

        for( int i = 0; i < nodesPerFace; i++) {
            un[2*i+0] = u[nodes[i]];
            un[2*i+1] = v[nodes[i]];
        }

        Matrix2d strain, stress;

        VectorXd vstrain;

        MatrixXd B;
        B.resize(3, dofPerFace);

        Point2D uv;
        for( int ig = 0; ig < numFaceGaussPoints; ig++) {
            uv[0] = fxi[ig];
            uv[1] = feta[ig];
            BMatrix(face, uv, B);
            vstrain = B*un;
            double exx = vstrain[0];
            double eyy = vstrain[1];
            double exy = 0.5*vstrain[2];
            double sxx = D(0,0)*exx + D(0,1)*eyy;
            double syy = D(1,0)*exx + D(1,1)*eyy;
            double sxy = 2.0*D(2,2)*exy;
            double maxVal = sqrt(sxx*sxx + syy*syy - sxx*syy + 3.0*sxy*sxy);
            if( maxVal > faceVonMises[fid] ) {
                faceStrain[fid] = strain;
                stress(0,0) = sxx;
                stress(0,1) = sxy;
                stress(1,0) = sxy;
                stress(1,1) = syy;
                faceStress[fid] = stress;
                faceVonMises[fid] = maxVal;
            }
        }
    */
}

//////////////////////////////////////////////////////////////////////////////////////////

void JElasticity2D:: computeStress()
{
    if(mesh == nullptr) return;

    if( verbose )
        cout << "Info: Computing face stress " << endl;

    size_t numnodes = mesh->getSize(0);
    size_t numfaces = mesh->getSize(2);

    nodeVonMises.resize(numnodes);
    for( size_t i = 0; i < numnodes; i++)
        nodeVonMises[i] = 0.0;

    faceVonMises.resize(numfaces);
    for( size_t i = 0; i < numfaces; i++)
        faceVonMises[i] = 0.0;

    faceStress.resize(numfaces);
    faceStrain.resize(numfaces);
    for( size_t i = 0; i < numfaces; i++) {
        faceStress[i].setZero();
        faceStrain[i].setZero();
    }

    // Calculate the strain and stress at every face ....
    for( size_t i = 0; i < numfaces; i++) {
        const JFacePtr &face = mesh->getFaceAt(i);
        computeStress(face);
    }

    if( verbose )
        cout << "Info: Computing nodes stress " << endl;


    // Calculate the Von-Mises stress on each nodes of the
    // mesh. ( Presently only on the origial nodes and not
    // on any higher order nodes ...)

    mesh->buildRelations(0,2);
    JFaceSequence faceneighs;
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vertex = mesh->getNodeAt(i);
        JVertex::getRelations(vertex, faceneighs);
        double sum = 0.0;
        int    numneighs = faceneighs.size();
        for( int j = 0; j < numneighs; j++) {
            int fid = faceneighs[j]->getID();
            sum += faceVonMises[fid];
        }
        nodeVonMises[i] = sum/(double)numneighs;
    }

    if( debugStep == 8) {
        cout << "VOnMises faces " << endl;
        latexVector( faceVonMises );
        cout << "VOnMises nodes " << endl;
        latexVector( nodeVonMises );
    }

    if( verbose)
        cout <<"Info: Stress calculations completed: All done " << endl;
}

/////////////////////////////////////////////////////////////////////////////////////
void JElasticity2D:: saveAs( const string &filename, const string &var)
{
    if( mesh == nullptr ) return;

    if( u.empty() || v.empty() ) return;

    size_t numnodes = mesh->getSize(0);
    JMeshVTKExporter mexp;

    if( var == "U" || var == "V") {
        if( !v.empty())  {
            for( size_t i = 0; i < numnodes; i++) {
                const JNodePtr &vtx = mesh->getNodeAt(i);
                if( var == "U")
                    vtx->setAttribute(var, u[i]);
                else
                    vtx->setAttribute(var, v[i]);
            }
        }
        mexp.addNodeAttribute(var);
    }

    if( var == "S") {
        computeStress();
        for( size_t i = 0; i < numnodes; i++) {
            const JNodePtr &vtx = mesh->getNodeAt(i);
            vtx->setAttribute(var, nodeVonMises[i]);
        }
        mexp.addNodeAttribute(var);
    }
    mexp.writeFile( mesh, filename);
}
#endif

/////////////////////////////////////////////////////////////////////////////////////

/*
void JElasticity2D :: concaveIntegrate(const JFacePtr &face, MatrixXd &Ke, vector<double> &fe)
{
    // Split the quad into two triangles: The concave corner becomes the divider...
    int corner = getConcaveCorner(face);
    if( corner < 0) return;

    JNodePtr p = face->getNodeAt(corner+2);
    face->setStartNode(p);

    Ke.setZero();
    fe.resize(dofPerFace);
    for( int i = 0; i < dofPerFace; i++)
        fe[i] = 0.0;

    concaveQuadTriIntegrate(face, 0, Ke, fe);
    concaveQuadTriIntegrate(face, 1, Ke, fe);

    if( debugStep == 2 ) {
        cout << "Calculating local stiffness matrix " << endl;
        latexMatrix(Ke);
        latexVector(fe);
    }
}
*/


