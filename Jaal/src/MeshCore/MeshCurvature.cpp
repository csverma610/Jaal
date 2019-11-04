#include "MeshCurvature.hpp"

#include <igl/gaussian_curvature.h>
#include <igl/principal_curvature.h>
#include <igl/jet.h>
#include <igl/avg_edge_length.h>

///////////////////////////////////////////////////////////////////////////////

void JMeshCurvature :: setMesh( const JMeshPtr &m)
{
    mesh = m;
    scale = 1.0;
    evalPos   = 0;

    if( mesh == nullptr) return;

    Array3D val;
    val[0] = 0.0;
    val[1] = 0.0;
    val[2] = 0.0;
    size_t numnodes = mesh->getSize(0);
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        vtx->setAttribute("Curvature", val);
    }
}
///////////////////////////////////////////////////////////////////////////////

int JMeshCurvature :: setGaussianCurvature()
{
    if( mesh == nullptr ) return 1;

    JMeshEigenMatrix mat;

    mat.setMesh(mesh);
    Eigen::MatrixXd V = mat.getNodeMatrix();
    Eigen::MatrixXi F = mat.getFaceMatrix();

    Eigen::VectorXd K;
    igl::gaussian_curvature(V,F,K);

    Array3D val;
    size_t numnodes = mesh->getSize(0);
    size_t index = 0;
    for( size_t i = 0; i < numnodes; i++) {
        const JNodePtr &vtx = mesh->getNodeAt(i);
        vtx->getAttribute("Curvature", K[index++]);
        val[0] = K[index++];
        vtx->setAttribute("Curvature", val);
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int JMeshCurvature :: setMeanCurvature( int method)
{
    if( mesh == nullptr) return 1;

    JMeshEigenMatrix mat;
    mat.setMesh(mesh);
    Eigen::MatrixXd V = mat.getNodeMatrix();
    Eigen::MatrixXi F = mat.getFaceMatrix();

    Eigen::VectorXd H;
    Eigen::MatrixXd PD1,PD2;
    Eigen::VectorXd PV1,PV2;
    Eigen::SparseMatrix<double> L,M,Minv;

    if( method == 0) {
        igl::principal_curvature(V,F,PD1,PD2,PV1,PV2);
        H = 0.5*(PV1+PV2);
    } else {
        igl::cotmatrix(V,F,L);
        igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_VORONOI,M);
        igl::invert_diag(M,Minv);
        // Laplace-Beltrami of position
        Eigen::MatrixXd HN = -Minv*(L*V);
        // Extract magnitude as mean curvature
        H = HN.rowwise().norm();
    }

    size_t numnodes = mesh->getSize(0);

    Array3D val;
    if( method == 0) {
        size_t index = 0;
        for( size_t i = 0; i < numnodes; i++) {
            const JNodePtr &vtx = mesh->getNodeAt(i);
            vtx->getAttribute("Curvature", val);
            val[1] = PV1[index];
            val[2] = PV2[index];
            vtx->setAttribute("Curvature", val);
            index++;
        }
    }

    if( evalPos == 2) {
        size_t numfaces = mesh->getSize(2);
        for( size_t i = 0; i < numfaces; i++) {
            const JFacePtr &f = mesh->getFaceAt(i);
            assignAverage(f, "Curvature");
        }
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
JMeshPtr JMeshCurvature :: getVectors( const Eigen::MatrixXd &headPoints,
                                       const Eigen::MatrixXd &tailPoints)
{
    JMeshPtr emesh = JMesh::newObject();
    Point3D  xyz;
    size_t numnodes = headPoints.rows();
    for( size_t i = 0; i < numnodes; i++) {
        JNodePtr  head = JNode::newObject();
        xyz[0] = headPoints.coeff(i,0);
        xyz[1] = headPoints.coeff(i,1);
        xyz[2] = headPoints.coeff(i,2);
        head->setXYZCoords( xyz );
        emesh->addObject( head );

        JNodePtr  tail = JNode::newObject();
        xyz[0] = tailPoints.coeff(i,0);
        xyz[1] = tailPoints.coeff(i,1);
        xyz[2] = tailPoints.coeff(i,2);
        tail->setXYZCoords( xyz );
        emesh->addObject( tail );

        JEdgePtr newedge = JEdge::newObject(head, tail);
        emesh->addObject(newedge);
    }
    return emesh;
}
///////////////////////////////////////////////////////////////////////////////
void JMeshCurvature :: assignAverage( const JFacePtr &face, const string &name)
{
    Array3D val, sum;

    sum[0] = 0.0;
    sum[1] = 0.0;
    sum[2] = 0.0;
    int   nn = face->getSize(0);
    for( int i = 0; i < nn; i++) {
        const JNodePtr &vtx = face->getNodeAt(i);
        vtx->getAttribute(name, val);
        sum[1] += val[1];
        sum[2] += val[2];
    }
    sum[1] /= (double)nn;
    sum[2] /= (double)nn;
    face->setAttribute(name, sum);
}
///////////////////////////////////////////////////////////////////////////////

std::pair<JMeshPtr, JMeshPtr> JMeshCurvature :: getCurvatureDirections()
{
    assert( mesh != nullptr);

    JMeshEigenMatrix mat;
    mat.setMesh(mesh);

    Eigen::MatrixXd V = mat.getNodeMatrix();
    Eigen::MatrixXi F = mat.getFaceMatrix();

    // Average edge length for sizing
    double avglen = 0.5*scale*igl::avg_edge_length(V,F);

    Eigen::MatrixXd PD1,PD2;
    Eigen::VectorXd PV1,PV2;

    igl::principal_curvature(V,F,PD1,PD2,PV1,PV2);

    Eigen::MatrixXd headPoints, tailPoints;

    headPoints = V + PD1*avglen;
    tailPoints = V - PD1*avglen;

    JMeshPtr minK = getVectors( headPoints, tailPoints);
    minK->setName("MinCurvatureVecField");

    headPoints = V + PD2*avglen;
    tailPoints = V - PD2*avglen;
    JMeshPtr maxK = getVectors( headPoints, tailPoints);
    maxK->setName("MaxCurvatureVecField");

    mesh->setAttribute("MinCurvatureVecField", minK);
    mesh->setAttribute("MaxCurvatureVecField", maxK);

    std::pair<JMeshPtr, JMeshPtr> minmax;
    minmax.first  = minK;
    minmax.second = maxK;
    return minmax;
}
///////////////////////////////////////////////////////////////////////////////

