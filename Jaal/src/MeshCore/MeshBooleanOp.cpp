#include "MeshBoolean.hpp"

/////////////////////////////////////////////////////////////////////////
int JMeshBoolean :: getOp( const string &str)
{
    if( str == "Union") return MESH_UNION;
    if( str == "Difference") return MESH_DIFFERENCE;
    if( str == "Intersection") return MESH_INTERSECTION;
    if( str == "Symmetric Difference") return MESH_SYMMETRIC_DIFFERENCE;
    return -1;
}

/////////////////////////////////////////////////////////////////////////

void JMeshBoolean :: setMesh( const JMeshPtr &meshA, const JMeshPtr &meshB)
{
    if( meshA == nullptr || meshB == nullptr)  return;

    JMeshEigenMatrix mat;

    mat.setMesh(meshA);
    VA = mat.getNodeMatrix();
    FA = mat.getFaceMatrix();

    mat.setMesh(meshB);
    VB = mat.getNodeMatrix();
    FB = mat.getFaceMatrix();
}

/////////////////////////////////////////////////////////////////////////


JMeshPtr  JMeshBoolean :: getMesh( int op )
{
    Eigen::VectorXi J;

    if( op == MESH_UNION) {
        igl::boolean::mesh_boolean(VA,FA,VB,FB,igl::boolean::MESH_BOOLEAN_TYPE_UNION,VC,FC, J);
    }

    if( op == MESH_INTERSECTION)
        igl::boolean::mesh_boolean(VA,FA,VB,FB,igl::boolean::MESH_BOOLEAN_TYPE_INTERSECT,VC,FC, J);

    if( op == MESH_DIFFERENCE)
        igl::boolean::mesh_boolean(VA,FA,VB,FB,igl::boolean::MESH_BOOLEAN_TYPE_MINUS,VC,FC, J);

    JMeshEigenMatrix mat;
    JMeshPtr meshC = mat.getMesh(VC, FC);

    int Aid = 0;
    int Bid = 1;
    for(size_t i = 0; i< FC.rows(); i++)
    {
        const JFacePtr &face = meshC->getFaceAt(i);
        if(J(i)<FA.rows())
            face->setAttribute("Partition", Aid);
        else
            face->setAttribute("Partition", Bid);
    }
    return meshC;
}
/////////////////////////////////////////////////////////////////////////
