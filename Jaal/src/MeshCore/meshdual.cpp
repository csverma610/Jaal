#include "MeshDUalGraph.hpp"

using namespace Jaal;

void usage()
{
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    Jaal::Mesh *mesh = new Jaal::Mesh2D;
    mesh->readFromFile( infilename );

    Jaal::MeshDualGraph mdual;
    mdual.setMesh(mesh);
    auto dualgraph = mdual.getGraph();

    mesh->saveAs( "dual.off");



    return 0;
}

