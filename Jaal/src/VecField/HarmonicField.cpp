#include "HarmonicField.hpp"

JHarmonicField :: JHarmonicField()
{
    igl::matlab::mlinit(&engine);
}

////////////////////////////////////////////////////////////////////////////////

JHarmonicField :: ~JHarmonicField()
{
    igl::matlab::mlclose(&engine);
}

////////////////////////////////////////////////////////////////////////////////

void JHarmonicField:: initSystem()
{
    if( mesh == nullptr) return;

    mesh->pruneNodes();
    mesh->pruneEdges();
    mesh->getTopology()->searchBoundary();

    mesh->getTopology()->getBoundary(boundNodes);
    mesh->getTopology()->getInternal(innerNodes);

    size_t numFree  = innerNodes.size();
    size_t numBound = boundNodes.size();
    if( numBound < 1 || numFree < 1) {
        cout << "Warning: There are no boundary or internal points " << endl;
        return;
    }

// enumerate nodes so that inner nodes are first;
    nodeLocalID.clear();
    size_t index = 0;
    for( const JNodePtr &vtx : innerNodes) nodeLocalID[vtx] = index++;
    for( const JNodePtr &vtx : boundNodes) nodeLocalID[vtx] = index++;

// Build Laplace Matrix :  Off Diagonal terms must be negative ,..
    size_t numNodes = mesh->getSize(0);
    size_t numEdges = mesh->getSize(1);
    JGeneralSparseMatrix<double> L;
    L.setSize(numNodes,numNodes);
    for( size_t i = 0; i < numEdges; i++) {
        const JEdgePtr &edge = mesh->getEdgeAt(i);
        assert( edge->isActive() );
        if( edge->isActive() ) {
            int v0 = nodeLocalID[edge->getNodeAt(0)];
            int v1 = nodeLocalID[edge->getNodeAt(1)];
            L.setValue(v0, v1, -1);
            L.setValue(v1, v0, -1);
        }
    }

// Diagonal terms must be positive ...
// JGeneralSparseMatrix<double>::SparseRow row;
    for( size_t i = 0; i < numNodes; i++) {
        const auto &row = L.getRow(i);
        assert( !row.empty() ) ;
        double sum = 0.0;
        for( auto  &keyVal: row)
            sum += keyVal.second;
        L.setValue(i,i, -sum);
    }

//  Now split the Matrix "L" into ( LA | Lb; 0 | I) parts...
    JGeneralSparseMatrix<double> Li, Lb;
    Li.setSize(numFree,numFree);
    Lb.setSize(numFree,numNodes-numFree);

    vector<size_t> rowIndex, colIndex,  scalarVal;
    for( size_t i = 0; i < numFree; i++) {
        const auto &row = L.getRow(i);
        for( auto  &keyVal: row) {
            size_t j   = keyVal.first;
            double val = keyVal.second;
            if( j < numFree )
                Li.setValue(i, j, val);
            else {
                Lb.setValue(i, j-numFree, val);
            }
        }
    }

    cout << "#Rows in the sparse matrix  :  " << Li.numRows() << endl;
    cout << "#Non-Zeros in sparse-matrix :  " << Li.getNumNonZeros() << endl;

//  Convert the Matrix Li, and Lb into Eigen Matrix,
    JEigenMatrixAdaptor  madpt;
    LA = madpt.getMatrix(Li, 1);
    LB = madpt.getMatrix(Lb, 1);

//  Set size of (xi,yi,zi) for internal nodes and (xb,yb,zb) for boundary nodes ...

//  Set the Matrix "LA" to the matlab, which will not change ....
    igl::matlab::mlsetmatrix(&engine,"A", LA);

//  Now initialize the solver to preprocess matrix "A"...
    switch( solver )
    {
    case CG:
        initCG();
        break;
    case CMG:
        initCMG();
        break;
    case HSC:
        initHSC();
        break;
    }
}

///////////////////////////////////////////////////////////////////////////////

int JHarmonicField:: solveSystem()
{
    JStopWatch  swatch;
    int status = -1;

    swatch.start();
    if( currStep == 0) initSystem();
    swatch.stop();
    double time1 = swatch.getSeconds();
    cout << "Harmonic field: Preprocessing time : " << time1  << endl;

    swatch.reset();

    swatch.start();
    switch(solver)
    {
    case CG:
        solveCG();
        break;
    case CMG:
        solveCMG();
        break;
    case HSC:
        solveHSC();
        break;
    }
    swatch.stop();

    double time2 = swatch.getSeconds();
    cout << "Harmonic field: calculations time  : " << time2  << endl;
    cout << "Total : " << time1 + time2 << endl;

    std::map<int,JNodePtr> nodeMap;
    for( const JNodePtr &vtx : innerNodes) {
        int id = nodeLocalID[vtx];
        nodeMap[id] = vtx;
    }

    size_t numFree = innerNodes.size();
    if( !scalarfield.empty() ) {
        for( size_t i = 0; i < numFree; i++) {
            double val = si.coeff(i,0);
            JNodePtr vtx = nodeMap[i];
            vtx->setAttribute(scalarfield, val);
        }
    } else {
        Point3D xyz;
        for( size_t i = 0; i < numFree; i++) {
            xyz[0] = xi.coeff(i,0);
            xyz[1] = yi.coeff(i,0);
            xyz[2] = zi.coeff(i,0);
            JNodePtr vtx = nodeMap[i];
            vtx->setXYZCoords(xyz);
        }
    }

    currStep++;
    return status;
}

///////////////////////////////////////////////////////////////////////////////

int JHarmonicField:: stdPCG( const string &cmd, Eigen::MatrixXd &xb, Eigen::MatrixXd &x)
{
    Eigen::MatrixXd b;
    int iter, flag, status = 0;
    double relres;

    // Solver for X-Coordinates ...
    b = -LB*xb;
    igl::matlab::mlsetmatrix(&engine,"b", b);
    igl::matlab::mleval(&engine, cmd);
    igl::matlab::mlgetmatrix(&engine,"X", x);
    flag   = igl::matlab::mlgetscalar(&engine, "flag");
    iter   = igl::matlab::mlgetscalar(&engine, "iter");
    relres = igl::matlab::mlgetscalar(&engine, "relres");
    status = max(status, flag);
    cout << "#Iterations : " << iter << " Final Residue " << relres << endl;
    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JHarmonicField:: stdPCG( const string &cmd)
{
    size_t numBound = boundNodes.size();
    size_t numFree  = innerNodes.size();

    if( !scalarfield.empty() )  {
        sb.resize(numBound,1);
        double val;
        for( size_t i = 0; i < numBound; i++) {
            boundNodes[i]->getAttribute(scalarfield, val);
            sb.coeffRef(i,0) = val;
        }
        si.resize(numFree,1);
        stdPCG(cmd, sb, si);
        return 0;
    }

    xb.resize(numBound,1);
    yb.resize(numBound,1);
    zb.resize(numBound,1);
    for( size_t i = 0; i < numBound; i++) {
        const Point3D &xyz = boundNodes[i]->getXYZCoords();
        xb.coeffRef(i,0) = xyz[0];
        yb.coeffRef(i,0) = xyz[1];
        zb.coeffRef(i,0) = xyz[2];
    }

    xi.resize(numFree,1);
    yi.resize(numFree,1);
    zi.resize(numFree,1);

    stdPCG(cmd, xb, xi);
    stdPCG(cmd, yb, yi);
    stdPCG(cmd, zb, zi);

    return 0;
}
///////////////////////////////////////////////////////////////////////////////

void JHarmonicField:: initCMG()
{
    // Combinatorial Multigrid method of Yiannis Koutis

    string cmd = "pfun = cmg_sdd(A)";
    igl::matlab::mleval(&engine, cmd);
}

///////////////////////////////////////////////////////////////////////////////
int JHarmonicField:: solveCMG()
{
    cout << "Solver : CMG " << endl;
    ostringstream oss;
    oss << "[X,flag,relres,iter] = pcg(A, b, " << tol << ", " << numIters << ", pfun)";
    string cmd = oss.str();
    return stdPCG(cmd);
}

///////////////////////////////////////////////////////////////////////////////
void JHarmonicField:: initCG()
{
//   string cmd = "L = ichol(A,struct('michol','on'))";
    string cmd = "L = ichol(A)";

    igl::matlab::mleval(&engine, cmd);
}

///////////////////////////////////////////////////////////////////////////////
int JHarmonicField:: solveCG()
{
    cout << "Solver : CG " << endl;
    ostringstream oss;
    oss << "[X,flag,relres,iter] = pcg(A, b, " << tol << "," << numIters << "," << "L, L')";
//  oss << "[X,flag,relres,iter] = pcg(A, b, " << tol << "," << numIters << ")";
    string cmd = oss.str();
    return stdPCG(cmd);
}
///////////////////////////////////////////////////////////////////////////////

void JHarmonicField:: initHSC()
{
    string cmd = "hsc_fun = hsc_setup(A,A)";
    igl::matlab::mleval(&engine, cmd);
}

///////////////////////////////////////////////////////////////////////////////
int JHarmonicField:: solveHSC()
{
    cout << "Solver : HSC " << endl;
    ostringstream oss;
    oss << "[X,flag,relres,iter] = pcg(A, b, " << tol << ", " << numIters << ", hsc_fun)";
    string cmd = oss.str();
    return stdPCG(cmd);
}

///////////////////////////////////////////////////////////////////////////////
