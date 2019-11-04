#include "MeshMeanCurvatureFlow.hpp"

JMeshMeanCurvatureFlow :: JMeshMeanCurvatureFlow()
{
    igl::matlab::mlinit(&engine);
}

JMeshMeanCurvatureFlow :: ~JMeshMeanCurvatureFlow()
{
    igl::matlab::mlclose(&engine);
}

void JMeshMeanCurvatureFlow :: setMesh(const JMeshPtr &m)
{
    orgMesh = m;
    currStep = 0;
    if( orgMesh == nullptr) return;

    int err = orgMesh->getAttribute("SmoothCurvatureMesh", mesh);
    if( err)  {
        mesh = orgMesh->deepCopy();
        orgMesh->setAttribute("SmoothCurvatureMesh", mesh);
        mesh->setName("Sphere");
    }

    affine.setMesh(mesh);
    affine.toCenter();
    maxDistance = 0.0;
    getSurfTriMesh();
    scaleFactor = 1.0;

}

///////////////////////////////////////////////////////////////////////////////

void JMeshMeanCurvatureFlow :: setAlgorithm( int a )
{
    algorithm = a;
    currStep  = 0;
}

///////////////////////////////////////////////////////////////////////////////

void JMeshMeanCurvatureFlow :: getSurfTriMesh()
{
    if( mesh == nullptr) return;

    int dim = mesh->getTopology()->getDimension();
    JMeshPtr  surfMesh;

    // If the mesh have 3 cells, extract the boundary mesh ...
    if( dim == 3)
        surfMesh = mesh->getTopology()->getSurfaceMesh();
    else
        surfMesh = mesh;

    orgBoundNodes = surfMesh->getNodes();

    triMesh.reset();
    AllTriMeshGenerator alltri;
    int elemType = surfMesh->getTopology()->isHomogeneous(2);
    if( elemType == JFace::QUADRILATERAL) {
        JMeshPtr tmpmesh = surfMesh->deepCopy();
        triMesh =  alltri.fromQuadMesh(tmpmesh,4);
    } else if( elemType == JFace::TRIANGLE)
        triMesh =  surfMesh;

    triMesh->enumerate(0);
    triBoundNodes = triMesh->getNodes();
    triMesh->getGeometry()->getCoordsArray( orgCoords, l2g);

    /*
        JDelaunayMesh2D del;
        if( intrinsicDelaunayMesh ) {
            del.setMesh(triMesh);
            del.getIntrinsicMesh();
        }
    */
}

////////////////////////////////////////////////////////////////////////////////

void JMeshMeanCurvatureFlow :: initMesh()
{
    if( triMesh == nullptr) return;

    surfArea[0] = triMesh->getGeometry()->getArea();
    sphRadius   = sqrt(surfArea[0]/(4.0*M_PI));

    triMesh->enumerate(0);
    triMesh->enumerate(1);
    triMesh->enumerate(2);

    currStep = 0;
    if( algorithm == IGL_METHOD) {
        JBoundingBox box = mesh->getGeometry()->getBoundingBox();
        scaleFactor  = box.getMaxLength();

        JMeshEigenMatrix mat;
        mat.setMesh(triMesh);
        V = mat.getNodeMatrix();
        F = mat.getFaceMatrix();

        mesh->enumerate(0);
        mesh->enumerate(2);

        getLMatrix();

        U = V;

        igl::doublearea(U,F,dblA);
        double area = 0.5*dblA.sum();
        igl::barycenter(U,F,BC);
        RowVector3d centroid(0,0,0);
        for(int i = 0; i<BC.rows(); i++)
        {
            centroid += 0.5*dblA(i)/area*BC.row(i);
        }
        U.rowwise() -= centroid;
        // Normalize to unit surface area (important for numerics)
        U.array() /= sqrt(area);
    }

    if( algorithm == KEENAN_METHOD)
        meshFair.setMesh(triMesh);
}
///////////////////////////////////////////////////////////////////////////////

void JMeshMeanCurvatureFlow ::restart()
{
    if( triMesh == nullptr) return;
    triMesh->enumerate(0);
    triMesh->getGeometry()->setCoordsArray(orgCoords, l2g);
    currStep = 0;

    /*
        if( algorithm == IGL_METHOD) {
            JMeshEigenMatrix mat;
            mat.setMesh(triMesh);
            U = mat.getNodeMatrix();
        }
    */
}
///////////////////////////////////////////////////////////////////////////////

void JMeshMeanCurvatureFlow :: updateMesh( const Eigen::MatrixXd &nodeCoords)
{
    size_t numnodes = orgBoundNodes.size();
    Point3D p0, p1;

    maxDistance = 0;
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = triBoundNodes[i];
        p0 = vtx->getXYZCoords();
        p1[0] = scaleFactor*nodeCoords.coeff(i,0);
        p1[1] = scaleFactor*nodeCoords.coeff(i,1);
        p1[2] = scaleFactor*nodeCoords.coeff(i,2);
        maxDistance = max(maxDistance, JMath::length2(p0,p1));
        orgBoundNodes[i]->setXYZCoords(p1);
    }
    maxDistance = sqrt( maxDistance );

    affine.toCenter();
    Point3D pc = mesh->getGeometry()->getCenter();
}

///////////////////////////////////////////////////////////////////////////////

void JMeshMeanCurvatureFlow :: predictTimeStep()
{
}

///////////////////////////////////////////////////////////////////////////////

int JMeshMeanCurvatureFlow:: stdPCG( const string &cmd)
{
    /*
        size_t numBound = boundNodes.size();
        for( size_t i = 0; i < numBound; i++) {
            const Point3D &xyz = boundNodes[i]->getXYZCoords();
            xb.coeffRef(i,0) = xyz[0];
            yb.coeffRef(i,0) = xyz[1];
            zb.coeffRef(i,0) = xyz[2];
        }

        Eigen::MatrixXd b;
        int iter, flag, status = 0;
        double relres;

        // Solver for X-Coordinates ...
        b = -LB*xb;
        igl::matlab::mlsetmatrix(&engine,"b", b);
        igl::matlab::mleval(&engine, cmd);
        igl::matlab::mlgetmatrix(&engine,"X", xi);
        flag   = igl::matlab::mlgetscalar(&engine, "flag");
        iter   = igl::matlab::mlgetscalar(&engine, "iter");
        relres = igl::matlab::mlgetscalar(&engine, "relres");
        status = max(status, flag);
        cout << "#Iterations : " << iter << " Final Residue " << relres << endl;

        // Solver for Y-Coordinates ...
        b = -LB*yb;
        igl::matlab::mlsetmatrix(&engine,"b", b);
        igl::matlab::mleval(&engine, cmd);
        igl::matlab::mlgetmatrix(&engine,"X", yi);
        iter   = igl::matlab::mlgetscalar(&engine, "iter");
        relres = igl::matlab::mlgetscalar(&engine, "relres");
        status = max(status, flag);
        cout << "#Iterations : " << iter << " Final Residue " << relres << endl;

        // Solver for Z-Coordinates ...
        b = -LB*zb;
        igl::matlab::mlsetmatrix(&engine,"b", b);
        igl::matlab::mleval(&engine, cmd);
        igl::matlab::mlgetmatrix(&engine,"X", zi);
        iter   = igl::matlab::mlgetscalar(&engine, "iter");
        relres = igl::matlab::mlgetscalar(&engine, "relres");
        status = max(status, flag);
        cout << "#Iterations : " << iter << " Final Residue " << relres << endl;
        return status;
    */

}
///////////////////////////////////////////////////////////////////////////////

void JMeshMeanCurvatureFlow :: solveSystem( const Eigen::SparseMatrix<double> &L,
        const Eigen::MatrixXd &b, Eigen::MatrixXd &x)
{
    if( linearSolver == CG_SOLVER) {
        Eigen::SimplicialLLT<Eigen::SparseMatrix<double > > solver(L);
        assert(solver.info() == Eigen::Success);
        x = solver.solve(b).eval();
        return;
    }


    /*
        JEigenMatrixAdpator madpt;
        JGeneralSparseMatrix<double>  A = madpt.getMatrix(L);
        int numRows = A.numRows();

        std::vector<Eigen::Triplet<double>>  triplets;
        triplets.reserve(Lt.getNumNonZeros);

        JGeneralSparseMatrix<int>::SparseRow row;
        double val;
        for( size_t i = 0; i < numRows; i++) {
            row = A.getRow(i);
            double sum = 0.0;
            for( auto &keyVal : row) {
                size_t j   = keyVal.first;
                val = keyVal.second;
                if( (i != j)  && (val < 0.0)) {
                    sum += val;
                    triplets.push_back(Eigen::Triplet<double>(i,j,val));
                }
            }
            assert( sum < 0.0);
            val = max(Lt.getVal(i,i), -sum);
            triplets.push_back(Eigen::Triplet<double>(i,i,val));
        }
        A.clear();
        Eigen::SparseMatrix<double> Lt;
        Lt.setFromTriplets( triplets.begin(), triplets.end() );

        igl::matlab::mlsetmatrix(&engine,"L", L);
        igl::matlab::mlsetmatrix(&engine,"Lt", Lt);

        string cmd = "hsc_fun = hsc_setup(Lt,L)";
        igl::matlab::mleval(&engine, cmd);
    */

    /*
        cout << "Solver : HSC " << endl;
        ostringstream oss;
        oss << "[X,flag,relres,iter] = pcg(A, b, " << tol << ", " << numIters << ", hsc_fun)";
        string cmd = oss.str();
        return stdPCG(cmd);
    */
}

///////////////////////////////////////////////////////////////////////////////
void JMeshMeanCurvatureFlow :: nextStep( int numIter)
{
    if( currStep == 0) initMesh();

    timeStep = max(1.0E-10, timeStep);
    if( algorithm == IGL_METHOD) {
        igl::massmatrix(U,F,igl::MASSMATRIX_TYPE_BARYCENTRIC,M);

        // Solve (M-delta*L) U = M*U
        const auto & S = (M - timeStep*L);
        Eigen::MatrixXd b = M*U;
        solveSystem(S, b, U);
        /*
                Eigen::SimplicialLLT<Eigen::SparseMatrix<double > > solver(S);
                assert(solver.info() == Eigen::Success);
                U = solver.solve(M*U).eval();
        */

        // Compute centroid and subtract (also important for numerics)
        igl::doublearea(U,F,dblA);
        double area = 0.5*dblA.sum();

        /*
                if( adaptiveTimeSteps ) {
                    if( fabs(area[1]-area[0])/fabs(area[0]) > 0.01) {
                        cout << "Reset time step " << endl;
                        timeStep /= 2.0;
                        return;
                    } else
                        timeStep *= 2.0;
                }
        */
        igl::barycenter(U,F,BC);
        centroid[0] = 0.0;
        centroid[1] = 0.0;
        centroid[2] = 0.0;

        for(int i = 0; i<BC.rows(); i++)
        {
            centroid += 0.5*dblA(i)/area*BC.row(i);
        }
        U.rowwise() -= centroid;
        // Normalize to unit surface area (important for numerics)
        U.array() /= sqrt(area);
        updateMesh(U);
    }

    if( algorithm == KEENAN_METHOD) {
        meshFair.setTimeStep(timeStep);
        meshFair.nextStep();
    }

    currStep++;
}

///////////////////////////////////////////////////////////////////////////////

void JMeshMeanCurvatureFlow :: projectOnSphere()
{
    if( mesh == nullptr) return;

    const Point3D &pc = mesh->getGeometry()->getCenter();

    size_t numnodes = mesh->getSize(0);

    Point3D pnew;
    Vec3D uvec;
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        const Point3D  &p1  = vtx->getXYZCoords();
        JMath::unit_vector(p1, pc, uvec);
        pnew[0] = pc[0] + sphRadius*uvec[0];
        pnew[1] = pc[1] + sphRadius*uvec[1];
        pnew[2] = pc[2] + sphRadius*uvec[2];
        vtx->setXYZCoords(pnew);
    }
}

///////////////////////////////////////////////////////////////////////////////

void JMeshMeanCurvatureFlow :: getLMatrix()
{
    Engine* engine;

    if(intrinsicDelaunayMesh) {
        igl::matlab::mlsetmatrix(&engine,"V",V);
        igl::matlab::mlsetmatrix(&engine,"F",F);
        ostringstream oss;
        oss << "[,F1] = intrinsic_delaunay_cotmatrix(V,F)";
        string cmd = oss.str();
        igl::matlab::mleval(&engine, cmd);
        igl::matlab::mlgetmatrix(&engine,"F1", F);
    }

    igl::cotmatrix(V,F,L);

    /*
       // Alternative construction of same Laplacian
        Eigen::SparseMatrix<double> G,K;
        Eigen::VectorXd dblA;
        // Gradient/Divergence
        igl::grad(V,F,G);
        // Diagonal per-triangle "mass matrix"
        igl::doublearea(V,F,dblA);
        // Place areas along diagonal #dim times
        const auto & T = 1.*(dblA.replicate(3,1)*0.5).asDiagonal();
        // Laplacian K built as discrete divergence of gradient or equivalently
        // discrete Dirichelet energy Hessian
        K = -G.transpose() * T * G;
    */

}


