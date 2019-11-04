#include "Mesh.hpp"

using namespace Jaal;

/////////////////////////////////////////////////////////////////////////////////

int main( int argc, char **argv)
{
/*
    JNodeSequence nodes(16);
    Mesh *msh = Mesh::newObject();

    for( int i = 0; i < 16; i++) {
        nodes[i] = Vertex::newObject();
        nodes[i]->setID(i);
    }
    msh->addObjects(nodes);
    JFaceSequence faces(9);

    int offset;
    for( int j = 0; j < 3; j++) {
        for( int i = 0; i < 3; i++) {
            offset = 4*j + i;
            Vertex *v0 = nodes[offset];
            Vertex *v1 = nodes[offset+1];
            offset = 4*(j+1) + i;
            Vertex *v2 = nodes[offset+ 1];
            Vertex *v3 = nodes[offset];
            Face *f = Quadrilateral::newObject(v0,v1,v2,v3);
            msh->addObject(f);
        }
    }

    cout << msh->getSize(0) << "  " << msh->getSize(1) << "  " << msh->getSize(2) << endl;
    cout << "*******************" << endl;

    for( int j = 0; j < 9; j++) {
        Face *f = msh->getFaceAt(j);
        f->setStatus( MeshEntity::REMOVE );
        cout << msh->getActiveSize(0) << "  " << msh->getActiveSize(1) << "  " << msh->getActiveSize(2) << endl;
        getchar();
    }

    exit(0);
*/

    int  xlength = 1.0;
    int  ylength = 1.0;
    int  npoints = 0;

    ofstream ofile( "test.poly", ios::out);

    ofile << "4 2 0 0" << endl;

    ofile << npoints   << "  " << -0.5*xlength << " " << -0.5*ylength << endl;
    ofile << npoints+1 << "  " <<  0.5*xlength << " " << -0.5*ylength << endl;
    ofile << npoints+2 << "  " <<  0.5*xlength << " " <<  0.5*ylength << endl;
    ofile << npoints+3 << "  " << -0.5*xlength << " " <<  0.5*ylength << endl;

    ofile << "4 1 " << endl;
    ofile << npoints   << "  " << npoints   << "  " << npoints+1 <<  " 1 " << endl;
    ofile << npoints+1 << "  " << npoints+1 << "  " << npoints+2 <<  " 2 " << endl;
    ofile << npoints+2 << "  " << npoints+2 << "  " << npoints+3 <<  " 3 " << endl;
    ofile << npoints+3 << "  " << npoints+3 << "  " << npoints   <<  " 4 " << endl;

    ofile << " 0 " << endl;
    ofile.close();

    system( "triangle -peq30a0.001 test.poly");

    Mesh *mesh = Mesh::newObject();
    mesh->readFromFile( "test.1.ele");
    mesh->saveAs("tri.xml");

    /*
       int dim[] = {10,5};
       Mesh *mesh = AllTriMeshGenerator::fromQuad4Tri( dim );
       mesh->saveAs("tri.vtk");
    */

    Mesh *quadmesh = AllQuadMeshGenerator::EdmondGraphMatching(mesh);
    cout << "Saving Quad mesh" << endl;
    quadmesh->saveAs("quad.xml");

    return 0;
    if( argc != 2) {
        cout << "Usage: executable <infile> <outfile> " << endl;
        return 1;
    }
}
